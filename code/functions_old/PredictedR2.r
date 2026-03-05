# PredictedR2.R: Cross-validated predicted R2 for mixed models
# Translation of MATLAB FREQUENTIST_ANALYSIS.m lines 501-585
library(data.table)
library(stats)
library(lme4)
library(lmerTest)
library(ggplot2)

# Clear workspace and graphics
rm(list = ls())
graphics.off()

# Source INITIALISE.R
if (!file.exists("code/INITIALISE.R")) stop("INITIALISE.R not found")
source("code/INITIALISE.R")

set.seed(1) # For reproducibility

# Outcomes that are positive indicators (higher = better wellbeing)
# These will be inverted to align with negative indicators
positive_outcomes <- c("SWL", "BRS", "FS")

# Define predictors (matching MixedModel.r)
all_predictors <- c("M_PosA", "M_NegA", "SD_PosA", "SD_NegA",
                    # "RMSSD_PosA", "RMSSD_NegA", "AR_PosA", "AR_NegA",
                    "P2N", "N2P")
bench <- c("P2N")
bench2plot <- "P2N"
ext <- c("P2N", "N2P")
ext2plot <- "P2N_N2P"

# Plot group names for consistent ordering
plot_group_names <- c("Well-being", "Depression", "Borderline", "Satisfaction")

# ============================================================================
# Helper Function: WithinBetweenSE
# Compute standard error of MSE using within and between group errors
# Based on formula 11.2 from Bayesian Data Analysis, Gelman 1992
# ============================================================================
WithinBetweenSE <- function(errors, groups) {
    uG <- unique(groups)
    n_groups <- length(uG)

    VW <- numeric(n_groups)     # within variance per group
    MSE_g <- numeric(n_groups)  # MSE per group
    W <- numeric(n_groups)      # group sizes

    for (i in seq_along(uG)) {
        idx <- which(groups == uG[i])
        VW[i] <- var(errors[idx])
        MSE_g[i] <- mean(errors[idx])
        W[i] <- length(idx)
    }

    MSE_weighted <- sum(W * MSE_g) / sum(W)
    Between <- sum(W * (MSE_g - MSE_weighted)^2) / sum(W)
    Within <- sum(W * VW) / sum(W)

    VarTotal <- Between / n_groups + Within / length(errors)
    return(sqrt(VarTotal))
}

# ============================================================================
# Helper Function: R2MixedCV
# 10-fold cross-validation for mixed models
# Returns squared prediction errors for each observation
# ============================================================================
R2MixedCV <- function(data, formula, group_var = "outcome_name_id", K = 10) {
    groups <- unique(data[[group_var]])
    n <- nrow(data)
    errors <- numeric(n)

    # Create fold assignments within each group (stratified CV)
    fold_idx <- integer(n)
    for (g in groups) {
        g_idx <- which(data[[group_var]] == g)
        n_g <- length(g_idx)
        # Assign folds cyclically within group
        fold_idx[g_idx] <- rep(1:K, length.out = n_g)
    }

    # Cross-validate across folds
    for (k in 1:K) {
        test_idx <- which(fold_idx == k)
        train_idx <- which(fold_idx != k)

        if (length(train_idx) < 10) next  # Skip if too few training samples

        # Fit mixed model on training data
        tryCatch({
            mdl <- lmer(formula, data = data[train_idx, ], REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 10000)))

            # Predict on test data (fixed effects only for out-of-sample prediction)
            pred <- predict(mdl, newdata = data[test_idx, ], re.form = NA,
                          allow.new.levels = TRUE)

            # Store squared errors
            errors[test_idx] <- (data$outcome_value[test_idx] - pred)^2
        }, error = function(e) {
            # On error, use mean prediction (baseline)
            errors[test_idx] <<- (data$outcome_value[test_idx] -
                                 mean(data$outcome_value[train_idx]))^2
        })
    }

    return(errors)
}

# ============================================================================
# Prepare data (same preprocessing as MixedModel.r)
# ============================================================================
dt <- long_data[outcome_type != "Classification", ]
outcome_types <- unique(dt$outcome_type)

cat("\n=== PREDICTED R2 ANALYSIS (CROSS-VALIDATED) ===\n")
cat(sprintf("Outcome types: %s\n", paste(outcome_types, collapse = ", ")))

dt[, outcome_type := as.factor(outcome_type)]
dt <- dt[!is.na(outcome_value)]
dt[outcome_name %in% positive_outcomes, outcome_value := -outcome_value]
dt[, outcome_name_id := rleid(outcome_name)]
dt[, outcome_value := scale(outcome_value)[, 1], by = .(dataset, outcome_name)]
dt[, (all_predictors) := lapply(.SD, function(x) scale(x)[, 1]),
   by = .(dataset, outcome_name), .SDcols = all_predictors]

# Filter to complete cases
initial_n <- nrow(dt)
complete_idx <- complete.cases(dt[, c("outcome_value", all_predictors), with = FALSE])
dt <- dt[complete_idx, ]
final_n <- nrow(dt)
cat(sprintf("Complete cases: %d -> %d rows (%.1f%% retained)\n\n",
           initial_n, final_n, 100 * final_n / initial_n))

# Create dummy variables for outcome_type
dummyOutcome <- model.matrix(~ 0 + outcome_type, data = dt)
outcome_type_levels <- levels(dt$outcome_type)

# Group variable for mixed models
GMixed <- dt$outcome_name_id
YMixed <- dt$outcome_value
var_Y <- var(YMixed)

# ============================================================================
# Results storage
# ============================================================================
results <- data.table(
    outcome_type = character(),
    model = character(),
    predictor = character(),
    R2_predicted = numeric(),
    SE = numeric()
)

# Total iterations for progress
total_preds <- length(all_predictors)
n_models <- 3  # solo, over_bench, over_ext
total_iter <- n_models * total_preds

cat("Starting cross-validated R2 analysis...\n")

# ============================================================================
# Loop over model types: solo (m=1), over_bench (m=2), over_ext (m=3)
# ============================================================================
iter <- 0
for (m in 1:3) {
    # Define basis for this model type
    if (m == 1) {
        basis <- c()
        model_name <- "solo"
    } else if (m == 2) {
        basis <- bench
        model_name <- "over_bench"
    } else {
        basis <- ext
        model_name <- "over_ext"
    }

    # Compute basis errors if basis is not empty
    if (length(basis) > 0) {
        # Well-being (no type interaction)
        formula_basis <- as.formula(paste0(
            "outcome_value ~ ", paste(basis, collapse = " + "),
            " + (1 + ", paste(basis, collapse = " + "), " | outcome_name_id)"
        ))
        errors_basis <- R2MixedCV(dt, formula_basis)

        # With type interaction (for per-type R2)
        # Build formula with dummy interactions
        fixed_terms <- c()
        for (b in basis) {
            for (j in seq_along(outcome_type_levels)) {
                fixed_terms <- c(fixed_terms,
                               paste0("I(", b, " * dummyOutcome[,", j, "])"))
            }
        }
        # Add intercept dummies
        for (j in seq_along(outcome_type_levels)) {
            fixed_terms <- c(fixed_terms, paste0("dummyOutcome[,", j, "]"))
        }
        formula_basis_type <- as.formula(paste0(
            "outcome_value ~ 0 + ", paste(fixed_terms, collapse = " + "),
            " + (1 + ", paste(basis, collapse = " + "), " | outcome_name_id)"
        ))
        errors_basis_type <- R2MixedCV(dt, formula_basis_type)
    }

    # Loop over predictors not in basis
    predictors_to_test <- setdiff(all_predictors, basis)

    for (pred in predictors_to_test) {
        iter <- iter + 1

        # Predictors for this model: basis + current predictor
        preds <- c(basis, pred)

        # ---- Well-being (no type interaction) ----
        formula_full <- as.formula(paste0(
            "outcome_value ~ ", paste(preds, collapse = " + "),
            " + (1 + ", paste(preds, collapse = " + "), " | outcome_name_id)"
        ))
        errors_full <- R2MixedCV(dt, formula_full)

        # Calculate R2 and SE
        if (length(basis) > 0) {
            # Added R2 over basis
            R2_pred <- (1 - mean(errors_full) / var_Y) -
                      (1 - mean(errors_basis) / var_Y)
            SE_pred <- WithinBetweenSE((errors_full - errors_basis)^2, GMixed) / var_Y
        } else {
            # Solo R2
            R2_pred <- 1 - mean(errors_full) / var_Y
            SE_pred <- WithinBetweenSE(errors_full, GMixed) / var_Y
        }

        results <- rbind(results, data.table(
            outcome_type = "Well-being",
            model = model_name,
            predictor = pred,
            R2_predicted = R2_pred,
            SE = SE_pred
        ))

        # ---- Per outcome type (with dummy interactions) ----
        # Build formula with dummy interactions for per-type R2
        fixed_terms <- c()
        for (p in preds) {
            for (j in seq_along(outcome_type_levels)) {
                fixed_terms <- c(fixed_terms,
                               paste0("I(", p, " * dummyOutcome[,", j, "])"))
            }
        }
        for (j in seq_along(outcome_type_levels)) {
            fixed_terms <- c(fixed_terms, paste0("dummyOutcome[,", j, "]"))
        }
        formula_full_type <- as.formula(paste0(
            "outcome_value ~ 0 + ", paste(fixed_terms, collapse = " + "),
            " + (1 + ", paste(preds, collapse = " + "), " | outcome_name_id)"
        ))
        errors_full_type <- R2MixedCV(dt, formula_full_type)

        # Calculate R2 for each outcome type
        for (t in seq_along(outcome_type_levels)) {
            otype <- outcome_type_levels[t]
            idx_t <- which(dummyOutcome[, t] == 1)
            var_Y_t <- var(YMixed[idx_t])

            if (length(basis) > 0) {
                R2_t <- (1 - mean(errors_full[idx_t]) / var_Y_t) -
                       (1 - mean(errors_basis[idx_t]) / var_Y_t)
                SE_t <- WithinBetweenSE(errors_full_type[idx_t] -
                                       errors_basis_type[idx_t], GMixed[idx_t]) / var_Y_t
            } else {
                R2_t <- 1 - mean(errors_full[idx_t]) / var_Y_t
                SE_t <- WithinBetweenSE(errors_full_type[idx_t], GMixed[idx_t]) / var_Y_t
            }

            results <- rbind(results, data.table(
                outcome_type = otype,
                model = model_name,
                predictor = pred,
                R2_predicted = R2_t,
                SE = SE_t
            ))
        }

        # Progress
        cat(sprintf("\r%.1f%% Completed (%s - %s)          ",
                   100 * iter / total_iter, model_name, pred))
    }
}

cat("\n\nAnalysis complete.\n")

# ============================================================================
# Save results
# ============================================================================
fwrite(results, "results/R2predicted.csv")
cat("Results saved to: results/R2predicted.csv\n")

# Print summary
cat("\n=== SUMMARY ===\n")
summary_dt <- results[, .(
    mean_R2 = mean(R2_predicted, na.rm = TRUE),
    sd_R2 = sd(R2_predicted, na.rm = TRUE),
    n = .N
), by = .(outcome_type, model)]
print(summary_dt)

# ============================================================================
# VISUALIZATION
# Following MATLAB FREQUENTIST_ANALYSIS.m lines 949-1002
# ============================================================================

# MATLAB colors
matlab_cols <- c("solo" = "#0072BD", "over_bench" = "#D95319", "over_ext" = "#EDB120")

# Create complete grid for consistent bar spacing
complete_grid <- CJ(
    predictor = all_predictors,
    model = c("solo", "over_bench", "over_ext"),
    outcome_type = plot_group_names
)

# Filter to valid model-predictor combinations
valid_combinations <- complete_grid[
    (model == "solo") |
        (model == "over_bench" & !(predictor %in% bench)) |
        (model == "over_ext" & !(predictor %in% ext))
]

# Left join to ensure all valid combinations exist
plot_data <- valid_combinations[results, on = .(predictor, model, outcome_type)]

# Make factors for correct order
plot_data[, predictor := factor(predictor, levels = all_predictors)]
plot_data[, model := factor(model, levels = c("solo", "over_bench", "over_ext"))]
plot_data[, outcome_type := factor(outcome_type, levels = plot_group_names)]

# Create faceted bar plot with error bars
p <- ggplot(plot_data, aes(x = predictor, y = R2_predicted, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7, preserve = "single"),
             width = 0.7, na.rm = TRUE) +
    geom_errorbar(aes(ymin = R2_predicted - SE, ymax = R2_predicted + SE),
                  width = 0.25, position = position_dodge(width = 0.7, preserve = "single"),
                  na.rm = TRUE) +
    scale_fill_manual(
        values = matlab_cols, name = "Model",
        labels = c("Solo", paste0("over ", bench2plot), paste0("over ", ext2plot))
    ) +
    labs(
        title = "Predicted R2 (Cross-Validated)",
        x = "Predictor",
        y = expression(Predicted ~ R^2)
    ) +
    facet_wrap(~outcome_type) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
    )

ggsave("plots/R2predicted.pdf", p, width = 12, height = 8)
cat("Plot saved to: plots/R2predicted.pdf\n")
