# Covariate analysis for solo models
library(data.table)
library(lme4)
library(lmerTest)

rm(list = ls())
graphics.off()

# Source INITIALISE.R (loads long_data, COVDataset, DataStrAll, etc.)
if (!file.exists("code/INITIALISE.R")) stop("INITIALISE.R not found")
source("code/INITIALISE.R")

set.seed(1)

# Outcomes that are positive indicators (higher = better wellbeing)
# These will be inverted to align with negative indicators (higher = worse wellbeing)
# for consistent interpretation of beta coefficients
positive_outcomes <- c("SWL", "BRS", "AAQ", "FS")

# --- PREPARE DATA ---
all_predictors <- c("M_PosA", "M_NegA", "SD_PosA", "SD_NegA",
                    "RMSSD_PosA", "RMSSD_NegA", "AR_PosA", "AR_NegA",
                    "P2N", "N2P")

dt <- long_data[outcome_type != "Classification", ]
dt[, outcome_type := as.factor(outcome_type)]
dt <- dt[!is.na(outcome_value)]
dt[outcome_name %in% positive_outcomes, outcome_value := -outcome_value] # invert positive indicators
dt[, outcome_name_id := rleid(outcome_name)]
dt[, outcome_value := scale(outcome_value)[, 1], by = .(dataset, outcome_name)]
dt[, (all_predictors) := lapply(.SD, function(x) scale(x)[, 1]), by = .(dataset, outcome_name), .SDcols = all_predictors]

# Merge dataset-level covariates to dt by "dataset"
covariate_names <- setdiff(names(wide_cov), "dataset")
dt <- merge(dt, wide_cov, by = "dataset", all.x = TRUE)

# Filter out rows with missing predictors or covariates once (more efficient than per-loop filtering)
complete_cols <- c(all_predictors, covariate_names, "outcome_value")
dt_complete <- dt[complete.cases(dt[, ..complete_cols])]

# --- COVARIATE ANALYSIS: SOLO MODELS ---
results_covariate <- data.table()

for (pred in all_predictors) {
    for (covariate in covariate_names) {
        # Fit the mixed-effects model
        # Use backticks to handle variable names that might contain special characters or reserved words
        # Include random slopes for predictor to match MATLAB implementation
        formula <- as.formula(
            paste0("outcome_value ~ `", pred, "` * `", covariate,
                   "` + (1 + `", pred, "` | outcome_name_id) + (1 | dataset)")
        )

        # Try to fit model with error handling
        mdl_result <- tryCatch({
            mdl <- lmer(formula, data = dt_complete)
            list(success = TRUE, model = mdl, error = NA, warning = NA, singular = isSingular(mdl))
        }, warning = function(w) {
            mdl <- lmer(formula, data = dt_complete)
            list(success = TRUE, model = mdl, error = NA, warning = w$message, singular = isSingular(mdl))
        }, error = function(e) {
            list(success = FALSE, model = NULL, error = e$message, warning = NA, singular = NA)
        })

        if (mdl_result$success) {
            s <- summary(mdl_result$model)
            coef_table <- s$coefficients

            # Find interaction term (handle both orderings: pred:cov or cov:pred)
            interaction_term <- NA_character_
            possible_terms <- c(paste0(pred, ":", covariate), paste0(covariate, ":", pred))
            for (term in possible_terms) {
                if (term %in% rownames(coef_table)) {
                    interaction_term <- term
                    break
                }
            }

            results_covariate <- rbind(
                results_covariate,
                data.table(
                    predictor = pred,
                    covariate = covariate,
                    beta = if (!is.na(interaction_term)) coef_table[interaction_term, "Estimate"] else NA_real_,
                    se = if (!is.na(interaction_term)) coef_table[interaction_term, "Std. Error"] else NA_real_,
                    t_value = if (!is.na(interaction_term)) coef_table[interaction_term, "t value"] else NA_real_,
                    p_value = if (!is.na(interaction_term)) coef_table[interaction_term, "Pr(>|t|)"] else NA_real_,
                    n_obs = nrow(mdl_result$model@frame),
                    converged = TRUE,
                    singular = mdl_result$singular,
                    warning = mdl_result$warning,
                    error = NA_character_
                )
            )
        } else {
            # Model failed to converge
            results_covariate <- rbind(
                results_covariate,
                data.table(
                    predictor = pred,
                    covariate = covariate,
                    beta = NA_real_,
                    se = NA_real_,
                    t_value = NA_real_,
                    p_value = NA_real_,
                    n_obs = NA_integer_,
                    converged = FALSE,
                    singular = NA,
                    warning = NA_character_,
                    error = mdl_result$error
                )
            )
        }
    }
}

# Apply multiple comparison corrections to p-values
# Only correct for converged models with non-missing p-values
converged_indices <- which(results_covariate$converged & !is.na(results_covariate$p_value))
n_tests <- length(converged_indices)

# Initialize correction columns with NA
results_covariate[, p_bonferroni := NA_real_]
results_covariate[, p_holm := NA_real_]
results_covariate[, p_fdr := NA_real_]

# Apply corrections only to converged models
if (n_tests > 0) {
    results_covariate[converged_indices, p_bonferroni := p.adjust(p_value, method = "bonferroni")]
    results_covariate[converged_indices, p_holm := p.adjust(p_value, method = "holm")]
    results_covariate[converged_indices, p_fdr := p.adjust(p_value, method = "fdr")]
}

# Save all results with descriptive statistics and corrected p-values
fwrite(results_covariate, file = "results/covariate_results.csv")

# Count significant interactions by different criteria
sig_uncorrected <- sum(results_covariate$p_value < 0.05 & results_covariate$converged == TRUE, na.rm = TRUE)
sig_bonferroni <- sum(results_covariate$p_bonferroni < 0.05 & results_covariate$converged == TRUE, na.rm = TRUE)
sig_holm <- sum(results_covariate$p_holm < 0.05 & results_covariate$converged == TRUE, na.rm = TRUE)
sig_fdr <- sum(results_covariate$p_fdr < 0.05 & results_covariate$converged == TRUE, na.rm = TRUE)

# Print summary
cat("\n=== COVARIATE ANALYSIS SUMMARY ===\n")
cat(sprintf("Total models fit: %d\n", nrow(results_covariate)))
cat(sprintf("Converged models: %d\n", sum(results_covariate$converged, na.rm = TRUE)))
cat(sprintf("Singular fits: %d\n", sum(results_covariate$singular, na.rm = TRUE)))
cat(sprintf("Models with warnings: %d\n", sum(!is.na(results_covariate$warning), na.rm = TRUE)))
cat(sprintf("Failed models: %d\n", sum(!results_covariate$converged, na.rm = TRUE)))
cat(sprintf("\nSignificant interactions (p < 0.05):\n"))
cat(sprintf("  - Uncorrected: %d\n", sig_uncorrected))
cat(sprintf("  - Bonferroni: %d\n", sig_bonferroni))
cat(sprintf("  - Holm: %d\n", sig_holm))
cat(sprintf("  - FDR (Benjamini-Hochberg): %d\n", sig_fdr))
cat(sprintf("\nResults saved to: results/covariate_results.csv\n"))