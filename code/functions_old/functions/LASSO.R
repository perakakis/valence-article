lassomodel_improved <- function(DV, predictors, nfolds_outer = 10, nfolds_inner = 10,
                                predictor_subsets = list(all = seq_len(ncol(predictors))),
                                use_loocv = FALSE, plot_results = TRUE) {
  # Input validation
  if (!is.numeric(DV)) stop("DV must be numeric")
  if (!is.matrix(predictors) && !is.data.frame(predictors)) stop("predictors must be a matrix or data frame")

  # Prepare data
  data <- data.frame(DV = DV, predictors)
  data <- na.omit(data)
  n <- nrow(data)
  y_var <- var(data$DV) # Variance of the dependent variable

  # Create outer folds
  if (use_loocv) {
    outer_folds <- as.list(seq_len(n)) # LOOCV: one fold per observation
    nfolds_outer <- n
  } else {
    outer_folds <- createFolds(data$DV, k = nfolds_outer, list = TRUE, returnTrain = FALSE)
  }

  # Initialize results
  results <- list()
  for (subset_name in names(predictor_subsets)) {
    subset_idx <- predictor_subsets[[subset_name]]
    cv_errors <- numeric(nfolds_outer)
    predicted_r2 <- numeric(nfolds_outer)
    selected_vars <- list()

    # Outer CV loop
    for (i in seq_along(outer_folds)) {
      test_idx <- outer_folds[[i]]
      train_data <- data[-test_idx, ]
      test_data <- data[test_idx, ]

      # Prepare training and test matrices
      x_train <- as.matrix(train_data[, subset_idx + 1]) # +1 to skip DV
      y_train <- train_data$DV
      x_test <- as.matrix(test_data[, subset_idx + 1])
      y_test <- test_data$DV

      # Inner CV for lambda selection
      cv_inner <- cv.glmnet(x_train, y_train,
        alpha = 1, nfolds = nfolds_inner,
        standardize = TRUE
      )

      # Fit model with selected lambda
      lasso_model <- glmnet(x_train, y_train,
        alpha = 1, lambda = cv_inner$lambda.min,
        standardize = TRUE
      )

      # Predict and compute error
      pred <- predict(lasso_model, newx = x_test, s = cv_inner$lambda.min)
      cv_errors[i] <- mean((y_test - pred)^2)

      # Compute predicted R-squared for this fold
      predicted_r2[i] <- 1 - cv_errors[i] / y_var

      # Track selected variables
      coefs <- coef(lasso_model)
      selected <- which(coefs[-1, ] != 0) # Exclude intercept
      selected_vars[[i]] <- selected
    }

    # Compute average CV error and predicted R-squared
    avg_cv_error <- mean(cv_errors)
    avg_predicted_r2 <- mean(predicted_r2)
    selection_freq <- sapply(seq_along(subset_idx), function(var) {
      mean(sapply(selected_vars, function(sel) var %in% sel))
    })

    # Fit final model on full data
    x_full <- as.matrix(data[, subset_idx + 1])
    y_full <- data$DV
    cv_full <- cv.glmnet(x_full, y_full,
      alpha = 1, nfolds = nfolds_inner,
      standardize = TRUE
    )
    final_model <- glmnet(x_full, y_full,
      alpha = 1, lambda = cv_full$lambda.min,
      standardize = TRUE
    )
    final_coefs <- coef(final_model)
    selected_predictors <- data.frame(
      predictor = colnames(predictors)[subset_idx],
      coefficient = as.numeric(final_coefs[-1, ])
    ) %>%
      filter(coefficient != 0)

    # Store results
    results[[subset_name]] <- list(
      final_model = final_model,
      selected_predictors = selected_predictors,
      cv_error = avg_cv_error,
      predicted_r2 = avg_predicted_r2,
      selection_frequency = data.frame(
        predictor = colnames(predictors)[subset_idx],
        frequency = selection_freq
      )
    )
  }

  # Plotting
  if (plot_results) {
    # Plot 1: Predicted R-squared across subsets
    r2_data <- data.frame(
      Subset = names(results),
      Predicted_R2 = sapply(results, function(x) x$predicted_r2)
    )
    p1 <- ggplot(r2_data, aes(x = Subset, y = Predicted_R2)) +
      geom_col(width = 0.2, fill = "#0072BD") +
      labs(
        title = "Predicted R-squared by Predictor Subset",
        x = "Subset", y = "Predicted R²"
      ) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5)
      )

    # Plot 2: Selection frequency for each subset
    plots_freq <- lapply(names(results), function(subset_name) {
      freq_data <- results[[subset_name]]$selection_frequency
      ggplot(freq_data, aes(x = predictor, y = frequency)) +
        geom_col(width = 0.3, fill = "#0072BD") +
        labs(
          title = paste("Variable Selection Frequency:", subset_name),
          x = "Predictor", y = "% Selected"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 16, hjust = 0.5)
        )
    })

    # Display plots
    print(p1)
    for (p in plots_freq) print(p)

    # Store plots in results
    results$plots <- list(predicted_r2 = p1, selection_frequency = plots_freq)
  }

  # Return results
  results
}