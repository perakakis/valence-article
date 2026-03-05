# RelativeImportance.R: Relative importance analysis for all continuous outcomes
library(relaimpo)   # For calculating relative importance
library(data.table)
library(ggplot2)
library(reshape2)

# Clear workspace and graphics
rm(list = ls())
graphics.off()

# Source INITIALISE.R
if (!file.exists("code/INITIALISE.R")) stop("INITIALISE.R not found")
source("code/INITIALISE.R")

set.seed(1) # For reproducibility

# Define predictors and outcome types (regression, not classification)
all_predictors <- c("M_PosA", "M_NegA", "SD_PosA", "SD_NegA",
                    "P2N", "N2P")
outcome_types <- c("Depression", "Anxiety", "Other")

# Print available outcomes by type
cat("\n=== RELATIVE IMPORTANCE ANALYSIS: OUTCOMES SUMMARY ===\n")
for (type in outcome_types) {
    outcomes_in_type <- unique(long_data[outcome_type == type, outcome_name])
    if (length(outcomes_in_type) > 0) {
        cat(sprintf("%s outcomes: %s\n", type, paste(outcomes_in_type, collapse = ", ")))
    }
}
cat("=======================================================\n\n")

# Embedded relative importance function (matches Spanish implementation)
relativeImportance <- function(DV, metrics) {
  data <- cbind(DV, metrics)
  data <- na.omit(data)

  full.model <- lm(DV ~ ., data = data)
  rel_imp <- calc.relimp(full.model, type = "lmg", rela = TRUE)
  return(rel_imp)
}

# Initialize results list
rel_imp_list <- list()

# Loop through datasets and outcomes
for (dataset_name in DataStrAll) {
    cat(sprintf("Processing dataset: %s\n", dataset_name))
    
    for (type in outcome_types) {
        outcomes_in_type <- unique(long_data[outcome_type == type, outcome_name])
        
        for (outcome in outcomes_in_type) {
            cat(sprintf("  Processing outcome: %s (%s)\n", outcome, type))
            
            dt <- long_data[dataset == dataset_name & outcome_name == outcome & outcome_type == type]
            if (nrow(dt) == 0) {
                cat("    No data found, skipping\n")
                next
            }

            # Aggregate to one row per participant 
            # Take first row per PID for outcome and all predictors
            cols_to_keep <- c("outcome_value", all_predictors)
            dt_agg <- dt[, lapply(.SD, function(x) head(x, 1)), by = PID, .SDcols = cols_to_keep]

            if (nrow(dt_agg) < 5) {
                cat("    Insufficient participants after aggregation, skipping\n")
                next
            }

            Y <- dt_agg$outcome_value
            metrics <- dt_agg[, ..all_predictors]

            # Perform relative importance analysis
            tryCatch({
                rel_imp <- relativeImportance(Y, metrics)

                # Calculate actual n_obs after na.omit inside relativeImportance
                data_for_n <- cbind(Y, metrics)
                n_obs <- nrow(na.omit(data_for_n))

                # Store results for each predictor
                rel_imp_list[[length(rel_imp_list) + 1]] <- data.table(
                    dataset = dataset_name,
                    outcome = outcome,
                    outcome_type = type,
                    predictor = names(rel_imp$lmg),
                    lmg_value = as.numeric(rel_imp$lmg),
                    n_obs = n_obs
                )
                cat(sprintf("    Completed: %d participants\n", n_obs))
            }, error = function(e) {
                cat(sprintf("    Error in relative importance analysis: %s\n", e$message))
            })
        }
    }
}

rel_imp_results <- rbindlist(rel_imp_list, use.names = TRUE)

# Check if we have any results
if (nrow(rel_imp_results) == 0) {
    stop("No relative importance results generated. Check data and outcome availability.")
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

# Calculate weighted means and SEs for relative importance by outcome_type and predictor
weighted_rel_imp <- rel_imp_results[, .(
    weighted_lmg = wmean(lmg_value, n_obs),
    weighted_lmg_se = wse(lmg_value, n_obs),
    n_total = sum(n_obs)
), by = .(outcome_type, predictor)]

# Create weighted results list for overall Well-being and per outcome_type
plot_group_names <- c("Well-being", outcome_types)
weighted_rel_imp_list <- list()

# Overall Well-being (aggregate across all outcome types)
weighted_rel_imp_list[["Well-being"]] <- rel_imp_results[, .(
    weighted_lmg = wmean(lmg_value, n_obs),
    weighted_lmg_se = wse(lmg_value, n_obs),
    n_total = sum(n_obs)
), by = .(predictor)]
weighted_rel_imp_list[["Well-being"]][, outcome_type := "Well-being"]

# Per outcome type
for (otype in outcome_types) {
    weighted_rel_imp_list[[otype]] <- weighted_rel_imp[outcome_type == otype, .(
        outcome_type, predictor, weighted_lmg, weighted_lmg_se, n_total
    )]
}

weighted_rel_imp_all <- rbindlist(weighted_rel_imp_list, use.names = TRUE, fill = TRUE)

# Save results
fwrite(rel_imp_results, "results/RelativeImportance.csv")
fwrite(weighted_rel_imp_all, "results/RelativeImportanceAvg.csv")

# MATLAB color (single color like LinearRegression.r)
matlab_color <- "#0072BD"

# Make factors for correct order
weighted_rel_imp_all[, predictor := factor(predictor, levels = all_predictors)]
weighted_rel_imp_all[, outcome_type := factor(outcome_type, levels = plot_group_names)]

# Create relative importance plot (similar to LinearRegression.r structure)
p <- ggplot(weighted_rel_imp_all, aes(x = predictor, y = weighted_lmg)) +
  geom_bar(stat = "identity", fill = matlab_color, width = 0.7) +
  geom_errorbar(aes(ymin = weighted_lmg - weighted_lmg_se, ymax = weighted_lmg + weighted_lmg_se),
                width = 0.25) +
  labs(
    title = "Relative Importance (LMG)",
    x = "Predictor",
    y = "Relative Importance"
  ) +
  facet_wrap(~outcome_type) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the plot
ggsave("plots/RelativeImportanceAvg.pdf", p, width = 12, height = 8)

# Create plot data for individual outcomes
individual_rel_imp_data <- rel_imp_results[, .(dataset, outcome, outcome_type, predictor, lmg_value, n_obs)]

# Calculate weighted relative importance for each individual outcome
individual_weighted_rel_imp <- individual_rel_imp_data[, .(
    weighted_lmg = wmean(lmg_value, n_obs),
    weighted_lmg_se = wse(lmg_value, n_obs),
    n_total = sum(n_obs)
), by = .(outcome, outcome_type, predictor)]

# Make factors for correct order
individual_weighted_rel_imp[, predictor := factor(predictor, levels = all_predictors)]

# Create individual outcomes plot
p_individual <- ggplot(individual_weighted_rel_imp, aes(x = predictor, y = weighted_lmg)) +
  geom_bar(stat = "identity", fill = matlab_color, width = 0.7) +
  geom_errorbar(aes(ymin = weighted_lmg - weighted_lmg_se, ymax = weighted_lmg + weighted_lmg_se),
                width = 0.25) +
  labs(
    title = "Relative Importance (LMG) by Individual Outcome",
    x = "Predictor",
    y = "Relative Importance"
  ) +
  facet_wrap(~outcome, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )
ggsave("plots/RelativeImportance.pdf", p_individual, width = 20, height = 12)

cat(sprintf("Generated %d relative importance results\n", nrow(rel_imp_results)))
cat(sprintf("Generated weighted averages for %d predictor-outcome combinations\n", nrow(weighted_rel_imp_all)))
cat("Results saved to:\n")
cat("  - results/RelativeImportance.csv\n")
cat("  - results/RelativeImportanceAvg.csv\n")
cat("  - plots/RelativeImportance.pdf\n")
cat("  - plots/RelativeImportanceAvg.pdf\n")