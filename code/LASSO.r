# LASSO (LOOCV, all predictors) for all continuous outcomes
library(data.table)
library(glmnet)
library(ggplot2)
library(R.utils)

# Clear workspace and graphics
rm(list = ls())
graphics.off()

# Define predictor sets for easy modification
all_predictors <- c("M_PosA", "M_NegA", "SD_PosA", "SD_NegA",
                    "P2N", "N2P")
bench_predictors <- c("P2N", "N2P")
ext_predictors <- c("P2N", "N2P", "SD_PosA", "SD_NegA")

# Define models for comparison
model_list <- list(
    all = all_predictors,
    bench = bench_predictors,
    ext = ext_predictors
)

# Path to cache file for results
cache_file <- "results/LassoCache.rds"
use_cache <- file.exists(cache_file)
outcome_types <- c("Depression", "Anxiety", "Other")
plot_group_names <- c("Well-being", outcome_types)

if (use_cache) {
    cache_data <- readRDS(cache_file)
    lasso_results <- cache_data$lasso_results
    select_freq <- cache_data$select_freq
    message("Loaded cached results from ", cache_file)
} else {
    # Source INITIALISE.R
    if (!file.exists("code/INITIALISE.R")) stop("INITIALISE.R not found")
    source("code/INITIALISE.R")

    set.seed(1) # For reproducibility

    # Store results
    lasso_results_list <- list()
    select_freq_list <- list()

    # Print available outcomes by type
    cat("\n=== LASSO ANALYSIS: OUTCOMES SUMMARY ===\n")
    for (type in outcome_types) {
        outcomes_in_type <- unique(long_data[outcome_type == type, outcome_name])
        if (length(outcomes_in_type) > 0) {
            cat(sprintf("%s outcomes: %s\n", type, paste(outcomes_in_type, collapse = ", ")))
        }
    }
    cat("=========================================\n\n")

    # Loop through datasets and outcomes
    total_datasets <- length(DataStrAll)
    dataset_count <- 0
    for (dataset_name in DataStrAll) {
        dataset_count <- dataset_count + 1
        cat(sprintf("Processing dataset %d/%d: %s\n", dataset_count, total_datasets, dataset_name))
        for (type in outcome_types) {
            outcomes_in_type <- unique(long_data[outcome_type == type, outcome_name])
            for (outcome in outcomes_in_type) {
                dt <- long_data[dataset == dataset_name & outcome_name == outcome & outcome_type == type]
                if (nrow(dt) == 0) next
                dt <- dt[!is.na(outcome_value)]
                if (nrow(dt) < 5) next # skip if too few

                Y <- scale(dt$outcome_value)[, 1]
                S <- as.data.table(scale(dt[, ..all_predictors]))

                # Filter to complete cases
                complete_idx <- complete.cases(S)
                n_complete <- sum(complete_idx)

                # Skip if too few complete cases
                if (n_complete < 5) {
                    cat(sprintf("  Skipping %s/%s: only %d complete cases\n",
                               dataset_name, outcome, n_complete))
                    next
                }

                Y <- Y[complete_idx]
                S <- S[complete_idx, ]
                n_obs <- length(Y)
                X <- as.matrix(S)

                # Compare all, bench, ext models for R2 (LOOCV for each set)
                for (model_name in names(model_list)) {
                    preds_r2 <- numeric(n_obs)
                    these_preds <- model_list[[model_name]]
                    X_model <- as.matrix(S[, ..these_preds])
                    for (i in seq_len(n_obs)) {
                        idx_train <- setdiff(seq_len(n_obs), i)
                        X_train <- X_model[idx_train, , drop = FALSE]
                        Y_train <- Y[idx_train]
                        X_test <- X_model[i, , drop = FALSE]
                        fit <- cv.glmnet(
                            X_train, Y_train,
                            alpha = 1, nfolds = min(20, length(Y_train)),
                            standardize = FALSE, intercept = TRUE
                        )
                        preds_r2[i] <- predict(fit, s = "lambda.min", newx = X_test)
                    }
                    ss_res <- sum((Y - preds_r2)^2)
                    ss_tot <- sum((Y - mean(Y))^2)
                    R2 <- 1 - ss_res / ss_tot

                    lasso_results_list[[length(lasso_results_list) + 1]] <- data.table(
                        dataset = dataset_name, outcome = outcome, outcome_type = type,
                        model = model_name, R2 = R2, n_obs = n_obs
                    )
                }

                # LOO-CV LASSO
                preds <- numeric(n_obs)
                selmat <- matrix(0, nrow = n_obs, ncol = length(all_predictors))
                colnames(selmat) <- all_predictors

                for (i in seq_len(n_obs)) {
                    idx_train <- setdiff(seq_len(n_obs), i)
                    X_train <- X[idx_train, , drop = FALSE]
                    Y_train <- Y[idx_train]
                    X_test <- X[i, , drop = FALSE]

                    # Fit LASSO with 20-fold CV (lambda selection inside fold)
                    fit <- cv.glmnet(
                        X_train, Y_train,
                        alpha = 1, nfolds = min(20, length(Y_train)),
                        standardize = FALSE, intercept = TRUE
                    )

                    # Predict left-out
                    preds[i] <- predict(fit, s = "lambda.min", newx = X_test)

                    # Which predictors were selected (nonzero beta at lambda.min)?
                    beta <- as.numeric(coef(fit, s = "lambda.min"))[-1] # drop intercept
                    selmat[i, ] <- as.integer(beta != 0)
                }

                # Selection frequency
                sel_freq <- colMeans(selmat)
                select_freq_list[[length(select_freq_list) + 1]] <- data.table(
                    dataset = dataset_name, outcome = outcome, outcome_type = type,
                    predictor = names(sel_freq), selection_freq = as.numeric(sel_freq),
                    n_obs = n_obs
                )
            }
        }
    }

    lasso_results <- rbindlist(lasso_results_list, use.names = TRUE)
    select_freq <- rbindlist(select_freq_list, use.names = TRUE)

    # Save cache for future runs
    saveRDS(list(lasso_results = lasso_results, select_freq = select_freq), cache_file)
    cat("Cache saved to: ", cache_file, "\n")
} # end of else

if (!exists("lasso_results") | !exists("select_freq")) {
    stop("Cache not loaded and results not available.")
}

# Weighted mean/SE functions
wmean <- function(x, w) sum(w * x, na.rm = TRUE) / sum(w[!is.na(x)])
wvar <- function(x, w) {
    wm <- wmean(x, w)
    w <- w[!is.na(x)]
    x <- x[!is.na(x)]
    sum_w <- sum(w)
    (sum(w * (x - wm)^2) / sum_w) * (sum_w / (sum_w - 1))
}
wse <- function(x, w) sqrt(wvar(x, w)) / sqrt(sum(!is.na(x)))

# R2: weighted averages by group and model
R2_weighted <- lasso_results[, .(
    weighted_R2 = wmean(R2, n_obs),
    weighted_R2_se = wse(R2, n_obs),
    n_total = sum(n_obs)
), by = .(outcome_type, model)]

# Well-being overall
wellbeing_R2 <- lasso_results[, .(
    weighted_R2 = wmean(R2, n_obs),
    weighted_R2_se = wse(R2, n_obs),
    n_total = sum(n_obs)
), by = model][, outcome_type := "Well-being"]

R2_weighted_all <- rbindlist(list(wellbeing_R2, R2_weighted), use.names = TRUE, fill = TRUE)
fwrite(R2_weighted_all, "results/LassoR2.csv")
cat("Results saved to: results/LassoR2.csv\n")

# MATLAB colors
matlab_blue <- "#0072BD"

# Plot grouped barplot: model comparison for each outcome_type
R2_weighted_all[, model := factor(model, levels = c("all", "bench", "ext"))]
R2_weighted_all[, outcome_type := factor(outcome_type, levels = plot_group_names)]
p_r2 <- ggplot(R2_weighted_all, aes(x = model, y = weighted_R2)) +
    geom_bar(
        stat = "identity", position = position_dodge(width = 0.7), width = 0.7,
        fill = matlab_blue
    ) +
    geom_errorbar(aes(ymin = weighted_R2 - weighted_R2_se, ymax = weighted_R2 + weighted_R2_se),
        width = 0.25, position = position_dodge(width = 0.7)
    ) +
    labs(
        title = "Weighted LASSO R² by Model and Outcome Type",
        x = "Model",
        y = expression(R^2)
    ) +
    facet_wrap(~outcome_type) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
ggsave("plots/LassoR2.pdf", p_r2, width = 12, height = 6)
cat("Plot saved to: plots/LassoR2.pdf\n")

# Selection frequency: weighted averages by group/predictor
select_weighted <- select_freq[, .(
    weighted_selection = wmean(selection_freq, n_obs),
    weighted_selection_se = wse(selection_freq, n_obs),
    n_total = sum(n_obs)
), by = .(outcome_type, predictor)]

# Well-being overall
wellbeing_select <- select_freq[, .(
    weighted_selection = wmean(selection_freq, n_obs),
    weighted_selection_se = wse(selection_freq, n_obs),
    n_total = sum(n_obs)
), by = predictor][, outcome_type := "Well-being"]

select_weighted_all <- rbindlist(list(wellbeing_select, select_weighted), use.names = TRUE, fill = TRUE)
fwrite(select_weighted_all, "results/LassoSelectionFreq.csv")
cat("Results saved to: results/LassoSelectionFreq.csv\n")

# Plot selection frequencies
select_weighted_all[, predictor := factor(predictor, levels = all_predictors)]
select_weighted_all[, outcome_type := factor(outcome_type, levels = plot_group_names)]
p_sel <- ggplot(select_weighted_all, aes(x = predictor, y = weighted_selection)) +
    geom_bar(
        stat = "identity", position = position_dodge(width = 0.7), width = 0.7,
        fill = matlab_blue, na.rm = TRUE
    ) +
    geom_errorbar(
        aes(
            ymin = weighted_selection - weighted_selection_se,
            ymax = weighted_selection + weighted_selection_se
        ),
        width = 0.25, position = position_dodge(width = 0.7), na.rm = TRUE
    ) +
    labs(
        title = "LASSO Selection Frequency by Outcome Type",
        x = "Predictor",
        y = "Selection Frequency"
    ) +
    facet_wrap(~outcome_type) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
ggsave("plots/LassoSelectionFreq.pdf", p_sel, width = 15, height = 8)
cat("Plot saved to: plots/LassoSelectionFreq.pdf\n")

# Create plots for individual outcomes (only outcomes present in data)
# R2 by individual outcome
individual_R2_weighted <- lasso_results[, .(
    weighted_R2 = wmean(R2, n_obs),
    weighted_R2_se = wse(R2, n_obs),
    n_total = sum(n_obs)
), by = .(outcome, outcome_type, model)]

individual_R2_weighted[, model := factor(model, levels = c("all", "bench", "ext"))]

p_r2_individual <- ggplot(individual_R2_weighted, aes(x = model, y = weighted_R2)) +
    geom_bar(
        stat = "identity", position = position_dodge(width = 0.7), width = 0.7,
        fill = matlab_blue
    ) +
    geom_errorbar(aes(ymin = weighted_R2 - weighted_R2_se, ymax = weighted_R2 + weighted_R2_se),
        width = 0.25, position = position_dodge(width = 0.7)
    ) +
    labs(
        title = "LASSO R² by Model and Individual Outcome",
        x = "Model",
        y = expression(R^2)
    ) +
    facet_wrap(~outcome, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10)
    )
ggsave("plots/LassoR2Individual.pdf", p_r2_individual, width = 20, height = 12)
cat("Plot saved to: plots/LassoR2Individual.pdf\n")

# Selection frequency by individual outcome
individual_select_weighted <- select_freq[, .(
    weighted_selection = wmean(selection_freq, n_obs),
    weighted_selection_se = wse(selection_freq, n_obs),
    n_total = sum(n_obs)
), by = .(outcome, outcome_type, predictor)]

individual_select_weighted[, predictor := factor(predictor, levels = all_predictors)]

fwrite(individual_select_weighted, "results/LassoSelectionFreqIndividual.csv")
cat("Results saved to: results/LassoSelectionFreqIndividual.csv\n")

# Predictor display labels with italic notation
predictor_labels <- c(
    "M_PosA" = "mPA", "M_NegA" = "mNA",
    "SD_PosA" = "sdPA", "SD_NegA" = "sdNA",
    "P2N" = "P2N", "N2P" = "N2P"
)
individual_select_weighted[, predictor_label := factor(
    predictor_labels[as.character(predictor)],
    levels = predictor_labels
)]

p_sel_individual <- ggplot(
    individual_select_weighted,
    aes(x = predictor_label, y = weighted_selection)
) +
    geom_col(
        width = 0.65, fill = "#4477AA", na.rm = TRUE
    ) +
    geom_errorbar(
        aes(
            ymin = pmax(weighted_selection - weighted_selection_se, 0),
            ymax = pmin(weighted_selection + weighted_selection_se, 1)
        ),
        width = 0.18, linewidth = 0.3, colour = "grey25", na.rm = TRUE
    ) +
    scale_y_continuous(
        limits = c(0, 1.08),
        breaks = seq(0, 1, 0.25),
        labels = c("0", ".25", ".50", ".75", "1"),
        expand = c(0, 0)
    ) +
    labs(
        x = NULL,
        y = "Selection Frequency"
    ) +
    facet_wrap(~outcome, ncol = 3) +
    theme_classic(base_size = 10) +
    theme(
        axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5, colour = "grey20"),
        axis.text.y = element_text(size = 7.5, colour = "grey20"),
        axis.title.y = element_text(size = 9.5, margin = margin(r = 6)),
        axis.line = element_line(linewidth = 0.3, colour = "grey40"),
        axis.ticks = element_line(linewidth = 0.25, colour = "grey40"),
        axis.ticks.length = unit(1.5, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(8, 10, 4, 4)
    )

n_outcomes <- length(unique(individual_select_weighted$outcome))
n_rows <- ceiling(n_outcomes / 3)
plot_height <- 2.0 + n_rows * 2.3
ggsave("plots/LassoSelectionFreqIndividual.pdf", p_sel_individual,
       width = 7, height = plot_height, dpi = 600)
cat("Plot saved to: plots/LassoSelectionFreqIndividual.pdf\n")