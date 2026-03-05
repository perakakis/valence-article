lassomodel <- function(DV, predictors) {
  # Input validation
  if (!is.numeric(DV)) stop("DV must be numeric")
  if (!is.matrix(predictors) && !is.data.frame(predictors)) stop("predictors must be a matrix or data frame")

  data <- data.frame(DV = DV, predictors)

# Omit rows with missing values
  initial_n <- nrow(data)
  data <- na.omit(data)
  final_n <- nrow(data)
  omitted_percentage <- (initial_n - final_n) / initial_n * 100
  cat("Omitted", round(omitted_percentage, 2), "% of data due to missing values.\n")

 
  y <- data$DV
  x <- model.matrix(DV ~ . - 1, data = data) # Design matrix without intercept

  # Cross-validated LASSO model
  cv.lasso <- cv.glmnet(x, y, alpha = 1, nfolds = 10, standardize = TRUE)

  # Fit final LASSO model using the lambda that minimized the cross-validation error
  lasso.model <- glmnet(x, y, alpha = 1, lambda = cv.lasso$lambda.min, standardize = TRUE)

  # Extract coefficients
  lasso.coefs <- coef(lasso.model)
  DV_lasso <- data.frame(
    predictor = row.names(lasso.coefs)[-1],
    coefficient = as.numeric(lasso.coefs[-1])
  )
  DV_lasso <- DV_lasso[DV_lasso$coefficient != 0, ]

  # Compute goodness-of-fit statistics
  predicted <- predict(lasso.model, newx = x, s = cv.lasso$lambda.min)
  r2 <- 1 - sum((y - predicted)^2) / sum((y - mean(y))^2)
  mse <- mean((y - predicted)^2)
  mae <- mean(abs(y - predicted))

  # Prepare result output
  result <- list(
    DV_lasso = DV_lasso,
    lambda_min = cv.lasso$lambda.min,
    cvm = min(cv.lasso$cvm), # Minimum mean cross-validated error
    cvsd = cv.lasso$cvsd[which.min(cv.lasso$cvm)], # Standard deviation of errors at min
    R2 = r2,
    MSE = mse,
    MAE = mae
  )

  return(result)
}
