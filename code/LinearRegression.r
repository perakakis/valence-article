# LinearRegression.R: Linear regression for all continuous outcomes
library(data.table)
library(stats)
library(caret)
library(ggplot2)
library(reshape2)

# Clear workspace and graphics
rm(list = ls())
graphics.off()

# Source INITIALISE.R
if (!file.exists("code/INITIALISE.R")) stop("INITIALISE.R not found")
source("code/INITIALISE.R")

set.seed(1) # For reproducibility

# Bootstrap configuration
n_boot <- 1000 # Number of bootstrap iterations
ci_alpha <- 0.05 # For 95% CI
cache_file <- "results/LinearRegressionCache.rds"

# Define predictors and outcome types
all_predictors <- c(
    "M_PosA", "M_NegA", "SD_PosA", "SD_NegA",
    "P2N", "N2P"
)

bench <- c("P2N")
bench2plot <- "P2N"
ext <- c("P2N", "N2P")
ext2plot <- "P2N_N2P"
outcome_types <- c("Depression", "Anxiety", "Other")

plot_group_names <- c("Well-being", outcome_types)

# Bootstrap helper function
bootstrap_R2 <- function(Y, S, bench, ext, all_predictors, n_boot) {
    n_obs <- length(Y)
    n_solo <- length(all_predictors)
    n_over_bench <- length(setdiff(all_predictors, bench))
    n_over_ext <- length(setdiff(all_predictors, ext))

    # Storage for bootstrap R2 values
    boot_solo <- matrix(NA, nrow = n_boot, ncol = n_solo)
    colnames(boot_solo) <- all_predictors
    boot_over_bench <- matrix(NA, nrow = n_boot, ncol = n_over_bench)
    colnames(boot_over_bench) <- setdiff(all_predictors, bench)
    boot_over_ext <- matrix(NA, nrow = n_boot, ncol = n_over_ext)
    colnames(boot_over_ext) <- setdiff(all_predictors, ext)

    for (b in 1:n_boot) {
        # Resample with replacement
        boot_idx <- sample(1:n_obs, n_obs, replace = TRUE)
        Y_boot <- Y[boot_idx]
        S_boot <- S[boot_idx, ]

        # Bench model R2 (needed for delta R2)
        X_bench <- as.matrix(S_boot[, bench, with = FALSE])
        mdl_bench <- lm(Y_boot ~ X_bench)
        R2_bench <- summary(mdl_bench)$r.squared

        # Ext model R2 (needed for delta R2)
        X_ext <- as.matrix(S_boot[, ext, with = FALSE])
        mdl_ext <- lm(Y_boot ~ X_ext)
        R2_ext <- summary(mdl_ext)$r.squared

        # Solo models
        for (i in seq_along(all_predictors)) {
            pred <- all_predictors[i]
            X <- S_boot[[pred]]
            mdl <- lm(Y_boot ~ X)
            boot_solo[b, i] <- summary(mdl)$r.squared
        }

        # Over bench models (delta R2)
        over_bench_preds <- setdiff(all_predictors, bench)
        for (i in seq_along(over_bench_preds)) {
            pred <- over_bench_preds[i]
            X <- as.matrix(S_boot[, c(bench, pred), with = FALSE])
            mdl <- lm(Y_boot ~ X)
            boot_over_bench[b, i] <- summary(mdl)$r.squared - R2_bench
        }

        # Over ext models (delta R2)
        over_ext_preds <- setdiff(all_predictors, ext)
        for (i in seq_along(over_ext_preds)) {
            pred <- over_ext_preds[i]
            X <- as.matrix(S_boot[, c(ext, pred), with = FALSE])
            mdl <- lm(Y_boot ~ X)
            boot_over_ext[b, i] <- summary(mdl)$r.squared - R2_ext
        }
    }

    list(solo = boot_solo, over_bench = boot_over_bench, over_ext = boot_over_ext)
}

# Print available outcomes by type
cat("\n=== LINEAR REGRESSION ANALYSIS: OUTCOMES SUMMARY ===\n")
for (type in outcome_types) {
    outcomes_in_type <- unique(long_data[outcome_type == type, outcome_name])
    if (length(outcomes_in_type) > 0) {
        cat(sprintf("%s outcomes: %s\n", type, paste(outcomes_in_type, collapse = ", ")))
    }
}
cat("=====================================================\n\n")

# Check for cached results
if (file.exists(cache_file)) {
    cat("Loading cached bootstrap results from:", cache_file, "\n")
    results <- readRDS(cache_file)
} else {
    # Initialize results list
    results_list <- list()

    # Loop through datasets and outcomes
    for (dataset_name in DataStrAll) {
        for (type in outcome_types) {
            outcomes_in_type <- unique(long_data[outcome_type == type, outcome_name])
            for (outcome in outcomes_in_type) {
                dt <- long_data[dataset == dataset_name & outcome_name == outcome & outcome_type == type]
                if (nrow(dt) == 0) next
                dt <- dt[!is.na(outcome_value)]
                Y <- scale(dt$outcome_value)[, 1]
                S <- as.data.table(scale(dt[, ..all_predictors]))

                # Filter to complete cases
                complete_idx <- complete.cases(S)
                n_complete <- sum(complete_idx)

                # Skip if too few complete cases
                if (n_complete < 5) {
                    cat(sprintf(
                        "  Skipping %s/%s: only %d complete cases\n",
                        dataset_name, outcome, n_complete
                    ))
                    next
                }

                Y <- Y[complete_idx]
                S <- S[complete_idx, ]
                n_obs <- length(Y)

                # Bootstrap confidence intervals
                cat(sprintf("  Bootstrapping %s/%s (%d obs)...\n", dataset_name, outcome, n_obs))
                boot_results <- bootstrap_R2(Y, S, bench, ext, all_predictors, n_boot)

                # Calculate percentile CIs for solo models
                solo_ci <- apply(boot_results$solo, 2, function(x) quantile(x, c(ci_alpha / 2, 1 - ci_alpha / 2)))

                # Calculate percentile CIs for over_bench models
                over_bench_ci <- apply(boot_results$over_bench, 2, function(x) quantile(x, c(ci_alpha / 2, 1 - ci_alpha / 2)))

                # Calculate percentile CIs for over_ext models
                over_ext_ci <- apply(boot_results$over_ext, 2, function(x) quantile(x, c(ci_alpha / 2, 1 - ci_alpha / 2)))

                # Bench model
                bench_str <- paste(bench, collapse = "_")
                X_bench <- as.matrix(S[, bench, with = FALSE])
                mdl_bench <- lm(Y ~ X_bench)
                summ_bench <- summary(mdl_bench)
                results_list[[length(results_list) + 1]] <- data.table(
                    dataset = dataset_name, outcome = outcome, outcome_type = type, model = "bench",
                    predictor = bench_str, R2 = summ_bench$r.squared, adj_R2 = summ_bench$adj.r.squared,
                    R2_over = 0, beta = NA, se = NA, p_value = NA,
                    n_obs = n_obs,
                    R2_ci_lower = NA, R2_ci_upper = NA,
                    R2_over_ci_lower = NA, R2_over_ci_upper = NA
                )

                # Ext model
                ext_str <- paste(ext, collapse = "_")
                X_ext <- as.matrix(S[, ext, with = FALSE])
                mdl_ext <- lm(Y ~ X_ext)
                summ_ext <- summary(mdl_ext)
                results_list[[length(results_list) + 1]] <- data.table(
                    dataset = dataset_name, outcome = outcome, outcome_type = type, model = "ext",
                    predictor = ext_str, R2 = summ_ext$r.squared, adj_R2 = summ_ext$adj.r.squared,
                    R2_over = 0, beta = NA, se = NA, p_value = NA,
                    n_obs = n_obs,
                    R2_ci_lower = NA, R2_ci_upper = NA,
                    R2_over_ci_lower = NA, R2_over_ci_upper = NA
                )

                # Solo models
                for (pred in all_predictors) {
                    X <- S[[pred]]
                    mdl <- lm(Y ~ X)
                    summ <- summary(mdl)
                    results_list[[length(results_list) + 1]] <- data.table(
                        dataset = dataset_name, outcome = outcome, outcome_type = type, model = "solo",
                        predictor = pred, R2 = summ$r.squared, adj_R2 = summ$adj.r.squared,
                        R2_over = 0, beta = coef(mdl)[2], se = coef(summ)[2, "Std. Error"],
                        p_value = coef(summ)[2, "Pr(>|t|)"],
                        n_obs = n_obs,
                        R2_ci_lower = solo_ci[1, pred], R2_ci_upper = solo_ci[2, pred],
                        R2_over_ci_lower = NA, R2_over_ci_upper = NA
                    )
                }

                # Over bench models
                for (pred in setdiff(all_predictors, bench)) {
                    X <- as.matrix(S[, c(bench, pred), with = FALSE])
                    mdl <- lm(Y ~ X)
                    summ <- summary(mdl)
                    idx <- length(bench) + 2
                    results_list[[length(results_list) + 1]] <- data.table(
                        dataset = dataset_name, outcome = outcome, outcome_type = type, model = "over_bench",
                        predictor = pred,
                        R2 = summ$r.squared, adj_R2 = summ$adj.r.squared,
                        R2_over = summ$r.squared - summ_bench$r.squared,
                        beta = coef(mdl)[idx], se = coef(summ)[idx, "Std. Error"],
                        p_value = coef(summ)[idx, "Pr(>|t|)"],
                        n_obs = n_obs,
                        R2_ci_lower = NA, R2_ci_upper = NA,
                        R2_over_ci_lower = over_bench_ci[1, pred], R2_over_ci_upper = over_bench_ci[2, pred]
                    )
                }

                # Over ext models
                for (pred in setdiff(all_predictors, ext)) {
                    X <- as.matrix(S[, c(ext, pred), with = FALSE])
                    mdl <- lm(Y ~ X)
                    summ <- summary(mdl)
                    idx <- length(ext) + 2
                    results_list[[length(results_list) + 1]] <- data.table(
                        dataset = dataset_name, outcome = outcome, outcome_type = type, model = "over_ext",
                        predictor = pred,
                        R2 = summ$r.squared, adj_R2 = summ$adj.r.squared,
                        R2_over = summ$r.squared - summ_ext$r.squared,
                        beta = coef(mdl)[idx], se = coef(summ)[idx, "Std. Error"],
                        p_value = coef(summ)[idx, "Pr(>|t|)"],
                        n_obs = n_obs,
                        R2_ci_lower = NA, R2_ci_upper = NA,
                        R2_over_ci_lower = over_ext_ci[1, pred], R2_over_ci_upper = over_ext_ci[2, pred]
                    )
                }
            }
        }
    }

    results <- rbindlist(results_list, use.names = TRUE)

    # Save to cache
    saveRDS(results, cache_file)
    cat("Bootstrap results cached to:", cache_file, "\n")
}

# Define weighted mean, variance and weighted SE functions
wmean <- function(x, w) sum(w * x, na.rm = TRUE) / sum(w[!is.na(x)])
wvar <- function(x, w) {
    wm <- wmean(x, w)
    w <- w[!is.na(x)]
    x <- x[!is.na(x)]
    sum_w <- sum(w)
    # Bessel's correction for weighted variance
    (sum(w * (x - wm)^2) / sum_w) * (sum_w / (sum_w - 1))
}
wse <- function(x, w) sqrt(wvar(x, w)) / sqrt(sum(!is.na(x)))

# Calculate weighted means and SEs for R2 and R2_over by outcome_type, predictor, and model
weighted_results <- results[, .(
    weighted_R2 = wmean(R2, n_obs),
    weighted_R2_se = wse(R2, n_obs),
    weighted_R2_ci_lower = wmean(R2_ci_lower, n_obs),
    weighted_R2_ci_upper = wmean(R2_ci_upper, n_obs),
    weighted_R2_over = wmean(R2_over, n_obs),
    weighted_R2_over_se = wse(R2_over, n_obs),
    weighted_R2_over_ci_lower = wmean(R2_over_ci_lower, n_obs),
    weighted_R2_over_ci_upper = wmean(R2_over_ci_upper, n_obs),
    n_total = sum(n_obs)
), by = .(outcome_type, predictor, model)]

# Create weighted_results_list for overall Well-being and per outcome_type
weighted_results_list <- list()
weighted_results_list[["Well-being"]] <- results[, .(
    weighted_R2 = wmean(R2, n_obs),
    weighted_R2_se = wse(R2, n_obs),
    weighted_R2_ci_lower = wmean(R2_ci_lower, n_obs),
    weighted_R2_ci_upper = wmean(R2_ci_upper, n_obs),
    weighted_R2_over = wmean(R2_over, n_obs),
    weighted_R2_over_se = wse(R2_over, n_obs),
    weighted_R2_over_ci_lower = wmean(R2_over_ci_lower, n_obs),
    weighted_R2_over_ci_upper = wmean(R2_over_ci_upper, n_obs),
    n_total = sum(n_obs)
), by = .(predictor, model)]

for (otype in outcome_types) {
    weighted_results_list[[otype]] <- weighted_results[outcome_type == otype, .(
        predictor, model,
        weighted_R2, weighted_R2_se,
        weighted_R2_ci_lower, weighted_R2_ci_upper,
        weighted_R2_over, weighted_R2_over_se,
        weighted_R2_over_ci_lower, weighted_R2_over_ci_upper,
        n_total
    )]
}

weighted_results_all <- rbindlist(weighted_results_list, idcol = "outcome_type", use.names = TRUE, fill = TRUE)

# Save weighted results
fwrite(weighted_results_all, "results/R2Avg.csv")
cat("Results saved to: results/R2Avg.csv\n")

# Extract R2 for each predictor (solo models) by dataset and outcome
solo_R2 <- results[model == "solo", .(dataset, outcome, predictor, R2)]
solo_R2 <- dcast(solo_R2, dataset + outcome ~ predictor, value.var = "R2")
solo_R2 <- as.data.table(solo_R2)
solo_R2[, ind := P2N > pmax(
    M_PosA, M_NegA, SD_PosA, SD_NegA, N2P,
    na.rm = TRUE
)]
rows <- solo_R2[ind == TRUE]
# print(rows)

fwrite(solo_R2, "results/R2.csv")
cat("Results saved to: results/R2.csv\n")

# MATLAB colors
matlab_cols <- c("solo" = "#0072BD", "over_bench" = "#D95319", "over_ext" = "#EDB120")

# Create complete grid to ensure consistent bar spacing
complete_grid <- CJ(
    predictor = all_predictors,
    model = c("solo", "over_bench", "over_ext"),
    outcome_type = plot_group_names
)

# Filter the complete grid to only include valid model-predictor combinations
valid_combinations <- complete_grid[
    (model == "solo") |
        (model == "over_bench" & !(predictor %in% bench)) |
        (model == "over_ext" & !(predictor %in% ext))
]

# Remove duplicate columns and filter out bench/ext models
weighted_results_clean <- weighted_results_all[
    !(model %in% c("bench", "ext")),
    !duplicated(names(weighted_results_all)),
    with = FALSE
]

# Left join to ensure all valid combinations exist (missing ones become NA)
plot_data <- valid_combinations[weighted_results_clean, on = .(predictor, model, outcome_type)]

# Make factors for correct order
plot_data[, predictor := factor(predictor, levels = all_predictors)]
plot_data[, model := factor(model, levels = c("solo", "over_bench", "over_ext"))]
plot_data[, outcome_type := factor(outcome_type, levels = plot_group_names)]

# Create plotting columns: use delta R2 for over models, raw R2 for solo
plot_data[, R2_plot := ifelse(
    model == "solo", weighted_R2, weighted_R2_over
)]
plot_data[, R2_se := ifelse(
    model == "solo", weighted_R2_se, weighted_R2_over_se
)]

# Create a single plot with consistent bar spacing
p <- ggplot(plot_data, aes(x = predictor, y = R2_plot, fill = model)) +
    geom_bar(
        stat = "identity", position = position_dodge(width = 0.7, preserve = "single"),
        width = 0.7, na.rm = FALSE
    ) +
    geom_errorbar(aes(ymin = R2_plot - R2_se, ymax = R2_plot + R2_se),
        width = 0.25, position = position_dodge(width = 0.7, preserve = "single"),
        na.rm = TRUE
    ) +
    scale_fill_manual(
        values = matlab_cols, name = "Model",
        labels = c("Solo", paste0("over ", bench2plot), paste0("over ", ext2plot))
    ) +
    labs(
        title = "Average R2",
        x = "Predictor",
        y = expression(Average ~ R^2)
    ) +
    facet_wrap(~outcome_type) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
    )
ggsave("plots/R2avg.pdf", p, width = 15, height = 8)
cat("Plot saved to: plots/R2avg.pdf\n")

# Create bootstrapped CI plot
# Use bootstrap percentile CIs: solo uses R2_ci, over models use R2_over_ci
plot_data[, R2_ci_lo := ifelse(
    model == "solo", weighted_R2_ci_lower, weighted_R2_over_ci_lower
)]
plot_data[, R2_ci_hi := ifelse(
    model == "solo", weighted_R2_ci_upper, weighted_R2_over_ci_upper
)]

p_ci <- ggplot(plot_data, aes(x = predictor, y = R2_plot, fill = model)) +
    geom_bar(
        stat = "identity", position = position_dodge(width = 0.7, preserve = "single"),
        width = 0.7, na.rm = FALSE
    ) +
    geom_errorbar(aes(ymin = R2_ci_lo, ymax = R2_ci_hi),
        width = 0.25, position = position_dodge(width = 0.7, preserve = "single"),
        na.rm = TRUE
    ) +
    scale_fill_manual(
        values = matlab_cols, name = "Model",
        labels = c("Solo", paste0("over ", bench2plot), paste0("over ", ext2plot))
    ) +
    labs(
        title = "Average R2 (Bootstrap 95% CI)",
        x = "Predictor",
        y = expression(Average ~ R^2)
    ) +
    facet_wrap(~outcome_type) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
    )
ggsave("plots/R2avg_CI.pdf", p_ci, width = 15, height = 8)
cat("Plot saved to: plots/R2avg_CI.pdf\n")

# Create plot data for individual outcomes
individual_plot_data <- results[
    !(model %in% c("bench", "ext")),
    .(dataset, outcome, outcome_type, model, predictor, R2, R2_over, n_obs)
]

# Calculate weighted R2 for each individual outcome
individual_weighted <- individual_plot_data[, .(
    weighted_R2 = wmean(R2, n_obs),
    weighted_R2_se = wse(R2, n_obs),
    weighted_R2_over = wmean(R2_over, n_obs),
    weighted_R2_over_se = wse(R2_over, n_obs),
    n_total = sum(n_obs)
), by = .(outcome, outcome_type, predictor, model)]

# Create complete grid for individual outcomes
individual_complete_grid <- CJ(
    predictor = all_predictors,
    model = c("solo", "over_bench", "over_ext"),
    outcome = unique(individual_weighted$outcome)
)

# Filter the complete grid to only include valid model-predictor combinations
individual_valid_combinations <- individual_complete_grid[
    (model == "solo") |
        (model == "over_bench" & !(predictor %in% bench)) |
        (model == "over_ext" & !(predictor %in% ext))
]

# Left join to ensure all valid combinations exist
individual_plot_final <- individual_valid_combinations[individual_weighted, on = .(predictor, model, outcome)]

# Make factors for correct order
individual_plot_final[, predictor := factor(predictor, levels = all_predictors)]
individual_plot_final[, model := factor(model, levels = c("solo", "over_bench", "over_ext"))]

# Create plotting columns: use delta R2 for over models, raw R2 for solo
individual_plot_final[, R2_plot := ifelse(
    model == "solo", weighted_R2, weighted_R2_over
)]
individual_plot_final[, R2_se := ifelse(
    model == "solo", weighted_R2_se, weighted_R2_over_se
)]

# Create individual outcomes plot
p_individual <- ggplot(individual_plot_final, aes(x = predictor, y = R2_plot, fill = model)) +
    geom_bar(
        stat = "identity", position = position_dodge(width = 0.7, preserve = "single"),
        width = 0.7, na.rm = FALSE
    ) +
    geom_errorbar(aes(ymin = R2_plot - R2_se, ymax = R2_plot + R2_se),
        width = 0.25, position = position_dodge(width = 0.7, preserve = "single"),
        na.rm = TRUE
    ) +
    scale_fill_manual(
        values = matlab_cols, name = "Model",
        labels = c("Solo", paste0("over ", bench2plot), paste0("over ", ext2plot))
    ) +
    labs(
        title = "Average R2 by Individual Outcome",
        x = "Predictor",
        y = expression(Average ~ R^2)
    ) +
    facet_wrap(~outcome, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        strip.text = element_text(size = 10)
    )
ggsave("plots/R2.pdf", p_individual, width = 20, height = 12)
cat("Plot saved to: plots/R2.pdf\n")


# P-VALUE AGGREGATION AND VISUALIZATION
# Calculate weighted p-values for ALL models by outcome_type
# Something to note: Weighted average of p-values is an approximation (Fisher's method would be more rigorous)

# Filter to models with p-values
p_results <- results[model %in% c("solo", "over_bench", "over_ext") & !is.na(p_value)]

# By outcome type
p_weighted_by_type <- p_results[, .(
    p_weighted = wmean(p_value, n_obs),
    n_total = sum(n_obs),
    n_outcomes = .N
), by = .(outcome_type, model, predictor)]

# Well-being
p_weighted_wellbeing <- p_results[, .(
    p_weighted = wmean(p_value, n_obs),
    n_total = sum(n_obs),
    n_outcomes = .N
), by = .(model, predictor)]
p_weighted_wellbeing[, outcome_type := "Well-being"]

# Combine
p_weighted_all <- rbind(
    p_weighted_wellbeing[, .(outcome_type, model, predictor, p_weighted, n_total, n_outcomes)],
    p_weighted_by_type[, .(outcome_type, model, predictor, p_weighted, n_total, n_outcomes)]
)

# Save p-value results
fwrite(p_weighted_all, "results/pValues.csv")
cat("P-values saved to: results/pValues.csv\n")

# Create p-value plot (log scale with significance thresholds)
cutoff <- 1e-6 # Floor for very small p-values
p_limit1 <- 0.05
p_limit2 <- 0.001

# Create complete grid for valid model-predictor combinations
p_complete_grid <- CJ(
    predictor = all_predictors,
    model = c("solo", "over_bench", "over_ext"),
    outcome_type = plot_group_names
)

# Filter to valid combinations
p_valid_combinations <- p_complete_grid[
    (model == "solo") |
        (model == "over_bench" & !(predictor %in% bench)) |
        (model == "over_ext" & !(predictor %in% ext))
]

# Merge with p-value data
p_plot_data <- p_valid_combinations[p_weighted_all, on = .(predictor, model, outcome_type)]
p_plot_data[, p_plot := pmax(p_weighted, cutoff)] # Apply floor
p_plot_data[, predictor := factor(predictor, levels = all_predictors)]
p_plot_data[, model := factor(model, levels = c("solo", "over_bench", "over_ext"))]
p_plot_data[, outcome_type := factor(outcome_type, levels = plot_group_names)]

# Create faceted p-value plot with all model types
p_pval <- ggplot(p_plot_data, aes(x = predictor, y = p_plot, fill = model)) +
    geom_bar(
        stat = "identity", position = position_dodge(width = 0.7, preserve = "single"),
        width = 0.7, na.rm = TRUE
    ) +
    geom_hline(yintercept = p_limit1, color = "#77AC30", linetype = "solid", linewidth = 0.8) +
    geom_hline(yintercept = p_limit2, color = "#77AC30", linetype = "dashed", linewidth = 0.8) +
    scale_fill_manual(
        values = matlab_cols, name = "Model",
        labels = c("Solo", paste0("over ", bench2plot), paste0("over ", ext2plot))
    ) +
    scale_y_log10(
        breaks = c(cutoff, p_limit2, p_limit1, 1),
        labels = c(
            paste0("<", format(cutoff, scientific = TRUE)),
            as.character(p_limit2), as.character(p_limit1), "1"
        )
    ) +
    labs(
        title = "P-values (Weighted Average)",
        x = "Predictor",
        y = "p-value"
    ) +
    facet_wrap(~outcome_type) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    )

ggsave("plots/pValues.pdf", p_pval, width = 12, height = 8)
cat("Plot saved to: plots/pValues.pdf\n")

# Also create individual outcome p-value plot
p_individual_outcomes <- results[model == "solo" & !is.na(p_value), .(
    p_weighted = wmean(p_value, n_obs),
    n_total = sum(n_obs)
), by = .(outcome, outcome_type, predictor)]

p_individual_outcomes[, p_plot := pmax(p_weighted, cutoff)]
p_individual_outcomes[, predictor := factor(predictor, levels = all_predictors)]

p_pval_individual <- ggplot(p_individual_outcomes, aes(x = predictor, y = p_plot)) +
    geom_bar(stat = "identity", fill = "#0072BD", width = 0.7) +
    geom_hline(yintercept = p_limit1, color = "#77AC30", linetype = "solid", linewidth = 0.8) +
    geom_hline(yintercept = p_limit2, color = "#77AC30", linetype = "dashed", linewidth = 0.8) +
    scale_y_log10(
        breaks = c(cutoff, p_limit2, p_limit1, 1),
        labels = c(
            paste0("<", format(cutoff, scientific = TRUE)),
            as.character(p_limit2), as.character(p_limit1), "1"
        )
    ) +
    labs(
        title = "Solo Variable P-values by Individual Outcome",
        x = "Predictor",
        y = "p-value"
    ) +
    facet_wrap(~outcome, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10)
    )

ggsave("plots/pValuesIndividual.pdf", p_pval_individual, width = 20, height = 12)
cat("Plot saved to: plots/pValuesIndividual.pdf\n")