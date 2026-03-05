library(data.table)
library(glmnet)
library(ggplot2)
library(pROC)
library(R.utils)

rm(list = ls())
graphics.off()

# Define predictor sets for easy modification
all_predictors <- c("M_PosA", "M_NegA", "SD_PosA", "SD_NegA",
                    "RMSSD_PosA", "RMSSD_NegA", "AR_PosA", "AR_NegA",
                    "P2N", "N2P")
bench_predictors <- c("P2N", "N2P")
ext_predictors <- c("P2N", "N2P", "SD_PosA", "SD_NegA")

# Outcomes to exclude (matching MATLAB LASSO.m line 231)
exclude_outcomes <- c("SCID_MDD_BPD")

# Define models for comparison
model_list <- list(
    all = all_predictors,
    bench = bench_predictors,
    ext = ext_predictors
)

cache_file <- "results/LassoClassCache.rds"
use_cache <- file.exists(cache_file)
outcome_types <- c("Classification")
plot_group_names <- c("Well-being", outcome_types)

if (use_cache) {
    cache_data <- readRDS(cache_file)
    lasso_results <- cache_data$lasso_results
    select_freq <- cache_data$select_freq
    message("Loaded cached results from ", cache_file)
} else {
    if (!file.exists("code/INITIALISE.R")) stop("INITIALISE.R not found")
    source("code/INITIALISE.R")

    set.seed(1)

    lasso_results <- data.table(
        dataset = character(), outcome = character(), model = character(), R2 = numeric(), n_obs = integer()
    )
    select_freq <- data.table()

    total_datasets <- length(DataStrAll)
    dataset_count <- 0
    for (dataset_name in DataStrAll) {
        dataset_count <- dataset_count + 1
        cat(sprintf("Processing dataset %d/%d: %s\n", dataset_count, total_datasets, dataset_name))
        outcomes_in_type <- unique(long_data[outcome_type == "Classification", outcome_name])
        # Exclude specified outcomes (matching MATLAB LASSO.m)
        outcomes_in_type <- setdiff(outcomes_in_type, exclude_outcomes)
        for (outcome in outcomes_in_type) {
            dt <- long_data[dataset == dataset_name & outcome_name == outcome & outcome_type == "Classification"]
            if (nrow(dt) == 0) next
            dt <- dt[!is.na(outcome_value)]
            Y <- dt$outcome_value
            # Check for valid binary class (0/1), and neither class is rare
            class_props <- prop.table(table(Y))
            if (!(length(class_props) == 2 && min(class_props) > 0.05)) next
            Y <- as.numeric(as.factor(Y)) - 1  # force 0/1
            S <- as.data.table(scale(dt[, ..all_predictors]))

            # Filter to complete cases
            complete_idx <- complete.cases(S)
            n_complete <- sum(complete_idx)

            # Skip if too few complete cases
            if (n_complete < 10) {
                cat(sprintf("  Skipping %s/%s: only %d complete cases\n",
                           dataset_name, outcome, n_complete))
                next
            }

            Y <- Y[complete_idx]
            S <- S[complete_idx, ]
            n_obs <- length(Y)

            for (model_name in names(model_list)) {
                these_preds <- model_list[[model_name]]
                X_model <- as.matrix(S[, ..these_preds])
                # Only for 'all' model: initialize selection matrix
                if (model_name == "all") {
                    selmat_all <- matrix(0, nrow = n_obs, ncol = length(these_preds))
                    colnames(selmat_all) <- these_preds
                }
                # LOO CV
                pred_prob <- numeric(n_obs)
                for (i in seq_len(n_obs)) {
                    idx_train <- setdiff(seq_len(n_obs), i)
                    X_train <- X_model[idx_train, , drop = FALSE]
                    Y_train <- Y[idx_train]
                    X_test <- X_model[i, , drop = FALSE]
                    fit <- glmnet(
                        X_train, Y_train,
                        family = "binomial",
                        alpha = 1,
                        standardize = FALSE,
                        intercept = TRUE
                    )
                    # Select lambda via cv.glmnet on training set (20-fold, matching MATLAB)
                    cvfit <- cv.glmnet(
                        X_train, Y_train,
                        family = "binomial",
                        alpha = 1,
                        nfolds = min(20, length(Y_train)),
                        standardize = FALSE,
                        intercept = TRUE
                    )
                    lambda_min <- cvfit$lambda.min
                    pred_prob[i] <- predict(fit, s = lambda_min, newx = X_test, type = "response")
                    if (model_name == "all") {
                        beta <- as.numeric(coef(fit, s = lambda_min))[-1]
                        selmat_all[i, ] <- as.integer(beta != 0)
                    }
                }
                # Compute accuracy-based R2 (matching MATLAB approach)
                # MATLAB: R2 = 1 - mean((round(pred_prob) - Y)^2)
                pred_class <- round(pred_prob)
                R2 <- 1 - mean((pred_class - Y)^2)
                lasso_results <- rbind(lasso_results, data.table(
                    dataset = dataset_name, outcome = outcome, model = model_name, R2 = R2, n_obs = n_obs
                ))
                if (model_name == "all") {
                    sel_freq <- colMeans(selmat_all)
                    select_freq <- rbind(
                        select_freq,
                        data.table(
                            dataset = dataset_name, outcome = outcome,
                            predictor = names(sel_freq), selection_freq = as.numeric(sel_freq),
                            n_obs = n_obs
                        )
                    )
                }
            }
        }
    }
    saveRDS(list(lasso_results = lasso_results, select_freq = select_freq), cache_file)
    cat("Cache saved to: ", cache_file, "\n")
}

if (!exists("lasso_results") | !exists("select_freq")) {
    stop("Cache not loaded and results not available.")
}

# Weighted mean, se
wmean <- function(x, w) sum(w * x, na.rm = TRUE) / sum(w[!is.na(x)])
wvar <- function(x, w) {
    wm <- wmean(x, w)
    w <- w[!is.na(x)]
    x <- x[!is.na(x)]
    sum_w <- sum(w)
    (sum(w * (x - wm)^2) / sum_w) * (sum_w / (sum_w - 1))
}
wse <- function(x, w) sqrt(wvar(x, w)) / sqrt(sum(!is.na(x)))

# Weighted R2 by model
r2_weighted <- lasso_results[, .(
    weighted_R2 = wmean(R2, n_obs),
    weighted_R2_se = wse(R2, n_obs),
    n_total = sum(n_obs)
), by = .(model)]

fwrite(r2_weighted, "results/LassoClassR2.csv")
cat("Results saved to: results/LassoClassR2.csv\n")

r2_weighted[, model := factor(model, levels = c("all", "bench", "ext"))]

# MATLAB colors
matlab_blue <- "#0072BD"

# Plot R2 comparison for all, bench, ext (classification)
r2_weighted[, model := factor(model, levels = c("all", "bench", "ext"))]
p_r2 <- ggplot(r2_weighted, aes(x = model, y = weighted_R2)) +
    geom_bar(stat = "identity", width = 0.7, fill = matlab_blue) +
    geom_errorbar(aes(ymin = weighted_R2 - weighted_R2_se, ymax = weighted_R2 + weighted_R2_se),
        width = 0.25) +
    labs(
        title = "Weighted LASSO R² by Model (Classification)",
        x = "Model",
        y = expression(R^2)
    ) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
ggsave("plots/LassoClassR2.pdf", p_r2, width = 8, height = 5)
cat("Plot saved to: plots/LassoClassR2.pdf\n")

if (nrow(select_freq) > 0) {
select_weighted <- select_freq[, .(
    weighted_selection = wmean(selection_freq, n_obs),
    weighted_selection_se = wse(selection_freq, n_obs),
    n_total = sum(n_obs)
), by = .(predictor)]

fwrite(select_weighted, "results/LassoClassSelectionFreq.csv")
cat("Results saved to: results/LassoClassSelectionFreq.csv\n")

# Plot selection frequencies
select_weighted[, predictor := factor(predictor, levels = all_predictors)]
p_sel <- ggplot(select_weighted, aes(x = predictor, y = weighted_selection)) +
    geom_bar(stat = "identity", width = 0.7, fill = matlab_blue, na.rm = TRUE) +
    geom_errorbar(aes(ymin = weighted_selection - weighted_selection_se, 
                     ymax = weighted_selection + weighted_selection_se),
        width = 0.25, na.rm = TRUE) +
    labs(
        title = "LASSO Variable Selection Frequency (Classification)",
        x = "Predictor",
        y = "Selection Frequency"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
ggsave("plots/LassoClassSelectionFreq.pdf", p_sel, width = 10, height = 5)
cat("Plot saved to: plots/LassoClassSelectionFreq.pdf\n")
} else {
    message("No selection frequency results to aggregate or plot.")
}