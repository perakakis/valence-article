# LinearRegression.R: Linear regression classification outcomes
library(data.table)
library(stats)
library(caret)
library(ggplot2)
library(pROC)

# Clear workspace and graphics
rm(list = ls())
graphics.off()

# Source INITIALISE.R
if (!file.exists("code/INITIALISE.R")) stop("INITIALISE.R not found")
source("code/INITIALISE.R")

set.seed(1) # For reproducibility

# Define predictors and outcome types (regression, not classification)
all_predictors <- c("M_PosA", "M_NegA", "SD_PosA", "SD_NegA",
                    "RMSSD_PosA", "RMSSD_NegA", "AR_PosA", "AR_NegA",
                    "P2N", "N2P")
# bench <- c("M_PosA", "M_NegA")
bench <- c("P2N")
bench2plot <- "P2N"
# ext <- c("M_PosA", "M_NegA", "SD_PosA", "SD_NegA")
ext <- c("P2N", "N2P")
ext2plot <- "P2N_N2P"
outcome_types <- c("Classification")

# plot_group_names <- c("Well-being", outcome_types)

# Initialize results data.table
results <- data.table(
    dataset = character(), outcome = character(), model = character(),
    predictor = character(), R2 = numeric(),
    R2_over = numeric(), beta = numeric(), se = numeric(), p_value = numeric(),
    n_obs = integer()
)

# Loop through datasets and outcomes
for (dataset_name in DataStrAll) {
    outcomes_in_type <- unique(long_data[outcome_type == "Classification", outcome_name])
    for (outcome in outcomes_in_type) {
        dt <- long_data[dataset == dataset_name & outcome_name == outcome & outcome_type == "Classification"]
        if (nrow(dt) == 0) next
        dt <- dt[!is.na(outcome_value)]
        Y <- dt$outcome_value
        # Check that both classes (0 and 1) are present and neither is too rare
        class_props <- prop.table(table(Y))
        if (!(length(class_props) == 2 && min(class_props) > 0.05)) next
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

        # Bench model
        bench_str <- paste(bench, collapse = "_")
        X_bench <- as.matrix(S[, bench, with = FALSE])
        mdl_bench <- glm(Y ~ X_bench, family = binomial(link = "logit"))
        auc_bench <- as.numeric(auc(Y, predict(mdl_bench, type = "response")))
        results <- rbind(results, data.table(
            dataset = dataset_name, outcome = outcome, model = "bench",
            predictor = bench_str, R2 = auc_bench,
            R2_over = 0, beta = NA, se = NA, p_value = NA,
            n_obs = n_obs
        ))

        # Ext model
        ext_str <- paste(ext, collapse = "_")
        X_ext <- as.matrix(S[, ext, with = FALSE])
        mdl_ext <- glm(Y ~ X_ext, family = binomial(link = "logit"))
        auc_ext <- as.numeric(auc(Y, predict(mdl_ext, type = "response")))
        results <- rbind(results, data.table(
            dataset = dataset_name, outcome = outcome, model = "ext",
            predictor = ext_str, R2 = auc_ext,
            R2_over = 0, beta = NA, se = NA, p_value = NA,
            n_obs = n_obs
        ))

        # Solo models
        for (pred in all_predictors) {
            X <- S[[pred]]
            mdl <- glm(Y ~ X, family = binomial(link = "logit"))
            auc <- as.numeric(auc(Y, predict(mdl, type = "response")))
            summ <- summary(mdl)
            results <- rbind(results, data.table(
                dataset = dataset_name, outcome = outcome, model = "solo",
                predictor = pred, R2 = auc,
                R2_over = 0, beta = coef(mdl)[2], se = coef(summ)[2, "Std. Error"],
                p_value = coef(summ)[2, "Pr(>|z|)"],
                n_obs = n_obs
            ))
        }

        # Over bench models
        for (pred in setdiff(all_predictors, bench)) {
            X <- as.matrix(S[, c(bench, pred), with = FALSE])
            mdl <- glm(Y ~ X, family = binomial(link = "logit"))
            auc <- as.numeric(auc(Y, predict(mdl, type = "response")))
            summ <- summary(mdl)
            idx <- length(bench) + 2
            results <- rbind(results, data.table(
                dataset = dataset_name, outcome = outcome, model = "over_bench",
                predictor = pred,
                R2 = auc,
                R2_over = auc - auc_bench,
                beta = coef(mdl)[idx], se = coef(summ)[idx, "Std. Error"],
                p_value = coef(summ)[idx, "Pr(>|z|)"],
                n_obs = n_obs
            ))
        }

        # Over ext models
        for (pred in setdiff(all_predictors, ext)) {
            X <- as.matrix(S[, c(ext, pred), with = FALSE])
            mdl <- glm(Y ~ X, family = binomial(link = "logit"))
            auc <- as.numeric(auc(Y, predict(mdl, type = "response")))
            summ <- summary(mdl)
            idx <- length(ext) + 2
            results <- rbind(results, data.table(
                dataset = dataset_name, outcome = outcome, model = "over_ext",
                predictor = pred,
                R2 = auc,
                R2_over = auc - auc_ext,
                beta = coef(mdl)[idx], se = coef(summ)[idx, "Std. Error"],
                p_value = coef(summ)[idx, "Pr(>|z|)"],
                n_obs = n_obs
            ))
        }
    }
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
    weighted_R2_over = wmean(R2_over, n_obs),
    weighted_R2_over_se = wse(R2_over, n_obs),
    n_total = sum(n_obs)
), by = .(predictor, model)]

# Save weighted results
fwrite(weighted_results, "results/R2ClassAvg.csv")
cat("Results saved to: results/R2ClassAvg.csv\n")

# Extract R2 for each predictor (solo models) by dataset and outcome
solo_R2 <- results[model == "solo", .(dataset, outcome, predictor, R2)]
solo_R2 <- dcast(solo_R2, dataset + outcome ~ predictor, value.var = "R2")
solo_R2 <- as.data.table(solo_R2)
solo_R2[, ind := P2N > pmax(
    M_PosA, M_NegA, SD_PosA, SD_NegA, N2P,
    na.rm = TRUE
)]
rows <- solo_R2[ind == TRUE]

fwrite(solo_R2, "results/R2Class.csv")
cat("Results saved to: results/R2Class.csv\n")

# MATLAB colors
matlab_cols <- c("solo" = "#0072BD", "over_bench" = "#D95319", "over_ext" = "#EDB120")

# Create complete grid to ensure consistent bar spacing
complete_grid <- CJ(
    predictor = all_predictors,
    model = c("solo", "over_bench", "over_ext")
)

# Filter the complete grid to only include valid model-predictor combinations
valid_combinations <- complete_grid[
    (model == "solo") |
        (model == "over_bench" & !(predictor %in% bench)) |
        (model == "over_ext" & !(predictor %in% ext))
]

# Remove duplicate columns and filter out bench/ext models 
weighted_results_clean <- weighted_results[
    !(model %in% c("bench", "ext")), 
    !duplicated(names(weighted_results)), 
    with = FALSE
]

# Left join to ensure all valid combinations exist (missing ones become NA)
plot_data <- valid_combinations[weighted_results_clean, on = .(predictor, model)]

# Make factors for correct order
plot_data[, predictor := factor(predictor, levels = all_predictors)]
plot_data[, model := factor(model, levels = c("solo", "over_bench", "over_ext"))]

# Create plotting columns: use delta R2 for over models, raw R2 for solo
plot_data[, R2_plot := ifelse(
    model == "solo", weighted_R2, weighted_R2_over
)]
plot_data[, R2_se := ifelse(
    model == "solo", weighted_R2_se, weighted_R2_over_se
)]

# Create a single plot with consistent bar spacing
p <- ggplot(plot_data, aes(x = predictor, y = R2_plot, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7, preserve = "single"), 
           width = 0.7, na.rm = FALSE) +
  geom_errorbar(aes(ymin = R2_plot - R2_se, ymax = R2_plot + R2_se),
                width = 0.25, position = position_dodge(width = 0.7, preserve = "single"), 
                na.rm = TRUE) +
  scale_fill_manual(
    values = matlab_cols, name = "Model",
    labels = c("Solo", paste0("over ", bench2plot), paste0("over ", ext2plot))
  ) +
  labs(
    title = "Average R2 for Classification Outcomes",
    x = "Predictor",
    y = expression(Average ~ R^2)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
ggsave("plots/R2ClassAvg.pdf", p, width = 15, height = 8)
cat("Plot saved to: plots/R2ClassAvg.pdf\n")

# Create plot data for individual outcomes
individual_plot_data <- results[
    !(model %in% c("bench", "ext")),
    .(dataset, outcome, model, predictor, R2, R2_over, n_obs)
]

# Calculate weighted R2 for each individual outcome
individual_weighted <- individual_plot_data[, .(
    weighted_R2 = wmean(R2, n_obs),
    weighted_R2_se = wse(R2, n_obs),
    weighted_R2_over = wmean(R2_over, n_obs),
    weighted_R2_over_se = wse(R2_over, n_obs),
    n_total = sum(n_obs)
), by = .(outcome, predictor, model)]

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
  geom_bar(stat = "identity", position = position_dodge(width = 0.7, preserve = "single"),
           width = 0.7, na.rm = FALSE) +
  geom_errorbar(aes(ymin = R2_plot - R2_se, ymax = R2_plot + R2_se),
                width = 0.25, position = position_dodge(width = 0.7, preserve = "single"),
                na.rm = TRUE) +
  scale_fill_manual(
    values = matlab_cols, name = "Model",
    labels = c("Solo", paste0("over ", bench2plot), paste0("over ", ext2plot))
  ) +
  labs(
    title = "Average R2 by Individual Classification Outcome",
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
ggsave("plots/R2Class.pdf", p_individual, width = 20, height = 12)
cat("Plot saved to: plots/R2Class.pdf\n")