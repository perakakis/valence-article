lassomodel <- function(DV, predictors) {
  data <- cbind(DV, predictors)

  initial_n <- nrow(data)
  data <- na.omit(data)
  final_n <- nrow(data)
  omitted_percentage <- (initial_n - final_n) / initial_n * 100
  cat("Omitted", round(omitted_percentage, 2), "% of data due to missing values.\n")

  y <- data$DV
  x <- model.matrix(DV ~ . - 1, data = data) # Remove intercept to fit LASSO
  x <- model.matrix(DV ~ ., data = data) # Remove intercept to fit LASSO

  # Cross-validated LASSO model
  cv.lasso <- cv.glmnet(x, y, alpha = 1, nfolds = 10)

  # Fit final LASSO model using the lambda that minimized the cross-validation error
  lasso.model <- glmnet(x, y, alpha = 1, lambda = cv.lasso$lambda.min)

  # Extract coefficients
  lasso.coefs <- coef(lasso.model)
  DV_lasso <- data.frame(
    predictor = row.names(lasso.coefs)[-1],
    coefficient = as.numeric(lasso.coefs[-1])
  )
  DV_lasso <- DV_lasso[DV_lasso$coefficient != 0, ]

  # Prepare result output
  result <- list(
    DV_lasso = DV_lasso,
    lambda_min = cv.lasso$lambda.min,
    cvm = min(cv.lasso$cvm), # Minimum mean cross-validated error
    cvsd = cv.lasso$cvsd[which.min(cv.lasso$cvm)] # Standard deviation of errors at min
  )

  return(result)
}

regressionmodel <- function(DV, predictor1, predictor2 = NULL) {
  if (!is.null(predictor2)) {
    data <- cbind(DV, predictor1, predictor2)
  } else {
    data <- cbind(DV, predictor1)
  }
  data <- data.frame(na.omit(data))

  # Create training control with 10-fold cross-validation
  ctrl <- trainControl(method = "cv", number = 10)

  # Fit linear regression with single predictor
  lm1 <- train(DV ~ predictor1,
    data = data, method = "lm",
    trControl = ctrl
  )

  # Save accuracy metrics
  allRsq1 <- lm1$resample$Rsquared
  Rsq1 <- lm1$results$Rsquared
  RMSE1 <- lm1$results$RMSE
  AIC1 <- AIC(lm1$finalModel)
  BIC1 <- BIC(lm1$finalModel)

  if (!is.null(predictor2)) {
    # Fit linear regression with two predictors
    lm2 <- train(DV ~ predictor1 + predictor2,
      data = data,
      method = "lm", trControl = ctrl
    )

    # Print accuracy metrics
    allRsq2 <- lm2$resample$Rsquared
    Rsq2 <- lm2$results$Rsquared
    RMSE2 <- lm2$results$RMSE
    AIC2 <- AIC(lm2$finalModel)
    BIC2 <- BIC(lm2$finalModel)

    # compare the performance of the two models using a likelihood ratio test (LRT)
    lrtest_res <- lrtest(lm1$finalModel, lm2$finalModel)
    pval <- lrtest_res$`Pr(>Chisq)`[2]

    return(list(
      allRsq1 = allRsq1, Rsq1 = Rsq1, RMSE1 = RMSE1, AIC1 = AIC1, BIC1 = BIC1,
      allRsq2 = allRsq2, Rsq2 = Rsq2, RMSE2 = RMSE2, AIC2 = AIC2, BIC2 = BIC2,
      pval = pval
    ))
  } else {
    return(list(allRsq1 = allRsq1, Rsq1 = Rsq1, RMSE1 = RMSE1, AIC1 = AIC1, BIC1 = BIC1))
  }
}

nestedmodels <- function(df, metrics, outcome) {
  library(lmtest)
  library(ggpubr)
  # Create an empty data frame to store the results
  results_df <- data.frame(
    result = numeric(length = length(metrics)),
    metric = metrics
  )

  # Loop over each metric for the current outcome
  for (i in 1:length(metrics)) {
    metric <- metrics[i]

    # Apply your function to the current metric and outcome
    result <- regressionmodel(df[[outcome]], df$P2N_ASR, df[[metric]])
    # result <- round(result$allRsq1,2)
    result <- round(result$pval, 3)

    # Add the result to the results data frame
    results_df[i, 1] <- result
  }

  my_theme <- list(fill = pal_jco("default")(2)[1], width = 0.5)

  ggbarplot(results_df,
    x = "metric", y = "result",
    color = pal_jco("default")(2)[1],
    fill = pal_jco("default")(2)[1],
    width = 0.3,
    ggpar = my_theme, position = position_dodge()
  ) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    ggtitle(outcome) +
    labs(x = "Metrics", y = "p-values")

  return(list(plot = last_plot(), results_over = results_df, outcome = outcome))
}

stepwiseRegression <- function(DV, metrics) {
  data <- cbind(DV, metrics)
  data <- na.omit(data)

  full.model <- lm(DV ~ ., data = data)
  step.model <- stepAIC(full.model,
    direction = "both",
    trace = FALSE
  )
  return(step.model)
}

relativeImportance <- function(DV, metrics) {
  data <- cbind(DV, metrics)

  initial_n <- nrow(data)
  data <- na.omit(data)
  final_n <- nrow(data)
  omitted_percentage <- (initial_n - final_n) / initial_n * 100
  cat("Omitted", round(omitted_percentage, 2), "% of data due to missing values.\n")

  full.model <- lm(DV ~ ., data = data)
  rel_imp <- calc.relimp(full.model,
    type = "lmg", rela = TRUE
  )
  # bootresults<-boot.relimp(full.model, b=1000)
  # ci<-booteval.relimp(bootresults, norank=T)
  # return(list(rel_imp = rel_imp, bootresults = bootresults,
  #             ci = ci))
  return(rel_imp)
}

identify_exclusions <- function(df, pid_col, val_col, duration_threshold) {
  # Initialize a vector to store the PIDs that meet the exclusion criteria
  exclusions <- c()

  # Unique participant IDs
  unique_pids <- unique(df[[pid_col]])

  # Loop over each unique participant
  for (pid in unique_pids) {
    # Subset data for the participant
    participant_data <- df[df[[pid_col]] == pid, ]

    # Initialize the flatline length counter
    flatline_length <- 0

    # Loop over the valence values
    for (i in seq_along(participant_data[[val_col]])) {
      # Check if the current value is -50, 0, 50, or NA
      if (is.na(participant_data[[val_col]][i]) || participant_data[[val_col]][i] %in% c(-50, 0, 50)) {
        # Increment the flatline counter
        flatline_length <- flatline_length + 1
      } else {
        # Reset the flatline counter if the current value is not -50, 0, 50, or NA
        flatline_length <- 0
      }

      # Check if the flatline length meets the threshold
      if (flatline_length >= duration_threshold) {
        exclusions <- c(exclusions, pid)
        break # Exit the loop for this participant since they met the exclusion criterion
      }
    }
  }

  # Return the list of participants to exclude
  return(unique(exclusions))
}

# Function to calculate reliability and produce comparison
calculate_reliability <- function(dfa, dfb, variable_name) {
  # Merge dataframes by 'PID'
  params_combined <- merge(dfa[, c("PID", variable_name)], dfb[, c("PID", variable_name)],
    by = "PID", suffixes = c(".half1", ".half2"), all = FALSE
  )

  # Correctly construct column names for filtering complete cases
  col_name_half1 <- paste(variable_name, ".half1", sep = "")
  col_name_half2 <- paste(variable_name, ".half2", sep = "")

  # Filter out rows where either half has NA for the given parameter
  complete_cases <- complete.cases(params_combined[, c(col_name_half1, col_name_half2)])
  params_filtered <- params_combined[complete_cases, ]

  # Calculate correlation and p-value if there are complete cases
  if (nrow(params_filtered) > 0) {
    test_results <- cor.test(params_filtered[[col_name_half1]], params_filtered[[col_name_half2]],
      method = "pearson", use = "complete.obs"
    )
    correlation <- test_results$estimate
    p_value <- test_results$p.value
  } else {
    correlation <- NA # Return NA if no complete cases
    p_value <- NA
  }

  # Optionally: Generate a plot to visualize the correlations
  if (nrow(params_filtered) > 0) {
    plot(params_filtered[[col_name_half1]], params_filtered[[col_name_half2]],
      main = paste("Correlation between halves for", variable_name),
      xlab = "First Half", ylab = "Second Half", pch = 19
    )
    abline(lm(params_filtered[[col_name_half2]] ~ params_filtered[[col_name_half1]]), col = "blue")
  }

  return(list(correlation = correlation, p_value = p_value))
}

# Function to calculate Mixed Emotions Index
calculate_mixed_emotions_index <- function(valence_data, range = c(-10, 10)) {
  valence_data <- valence_data[!is.na(valence_data)] # Remove NA values
  within_range <- valence_data >= range[1] & valence_data <= range[2]
  return(sum(within_range) / length(valence_data))
}
