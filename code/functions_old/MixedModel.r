# MixedModel.R: Mixed model analysis for continuous outcomes
library(data.table)
library(stats)
library(ggplot2)
library(lme4)
library(lmerTest)

# Clear workspace and graphics
rm(list = ls())
graphics.off()

# Source INITIALISE.R
if (!file.exists("code/INITIALISE.R")) stop("INITIALISE.R not found")
source("code/INITIALISE.R")

set.seed(1) # For reproducibility

# Outcomes that are positive indicators (higher = better wellbeing)
# These will be inverted to align with negative indicators (higher = worse wellbeing)
# for consistent interpretation of beta coefficients
positive_outcomes <- c("SWL", "BRS", "FS")

# Define predictors and outcome types
all_predictors <- c("M_PosA", "M_NegA", "SD_PosA", "SD_NegA",
                    "RMSSD_PosA", "RMSSD_NegA", "AR_PosA", "AR_NegA",
                    "P2N", "N2P")
bench <- c("P2N")
bench2plot <- "P2N"
ext <- c("P2N", "N2P")
ext2plot <- "P2N_N2P"

# Aux function to calculate R2 for mixed models
R2Mixed <- function(mod, types = NULL) {
    # Fitted values from fixed effects only
    y_fixed <- predict(mod, re.form = NA)
    # Residuals (includes both fixed and random)
    res <- residuals(mod)
    # Fitted values from full model
    y_fitted <- fitted(mod)
    # Residual variance (model error)
    se_vec <- res^2
    # Fixed effects variance (centered)
    sf_vec <- (y_fixed - mean(y_fixed))^2
    # Random effects variance (difference between fitted and fixed)
    sl_vec <- (y_fitted - y_fixed)^2
    # If no types, just overall R2
    if (is.null(types)) {
        se <- mean(se_vec)
        sf <- mean(sf_vec)
        sl <- mean(sl_vec)
        R2 <- (sl + sf) / (sl + sf + se)
        return(list(R2_mixed = R2))
    }
    # If types, calculate for each unique type
    type_levels <- unique(types)
    R2_type <- numeric(length(type_levels) + 1)
    # First, overall R2
    se <- mean(se_vec)
    sf <- mean(sf_vec)
    sl <- mean(sl_vec)
    R2_type[1] <- (sl + sf) / (sl + sf + se)
    # Then, for each type
    for (t in seq_along(type_levels)) {
        idx <- which(types == type_levels[t])
        se <- mean(se_vec[idx])
        sf <- mean(sf_vec[idx])
        sl <- mean(sl_vec[idx])
        R2_type[t + 1] <- (sl + sf) / (sl + sf + se)
    }
    names(R2_type) <- c("Well-being", as.character(type_levels))
    return(R2_type)
}

# Prepare data table
dt <- long_data[outcome_type != "Classification", ] # remove classification outcome
outcome_types <- c(unique(dt$outcome_type))

# Print available outcomes by type
cat("\n=== MIXED MODEL ANALYSIS: OUTCOMES SUMMARY ===\n")
for (type in outcome_types) {
    outcomes_in_type <- unique(dt[outcome_type == type, outcome_name])
    if (length(outcomes_in_type) > 0) {
        cat(sprintf("%s outcomes: %s\n", type, paste(outcomes_in_type, collapse = ", ")))
    }
}
cat("===============================================\n\n")

dt[, outcome_type := as.factor(outcome_type)]
dt <- dt[!is.na(outcome_value)] # remove NA
dt[outcome_name %in% positive_outcomes, outcome_value := -outcome_value] # invert positive indicators
dt[, outcome_name_id := rleid(outcome_name)] # create outcome_name_id
dt[, outcome_value := scale(outcome_value)[, 1], by = .(dataset, outcome_name)] # scale outcome values
dt[, (all_predictors) := lapply(.SD, function(x) scale(x)[, 1]),
    by = .(dataset, outcome_name), .SDcols = all_predictors
] # scale predictors

# Filter to complete cases
initial_n <- nrow(dt)
complete_idx <- complete.cases(dt[, c("outcome_value", all_predictors), with = FALSE])
dt <- dt[complete_idx, ]
final_n <- nrow(dt)
cat(sprintf("Complete cases filtering: %d -> %d rows (%.1f%% retained)\n",
           initial_n, final_n, 100 * final_n / initial_n))

# Results data table
mixed_results <- data.table(
    outcome_type = character(),
    model = character(),
    predictor = character(),
    R2 = numeric(),
    R2_over = numeric(),
    beta = numeric(),
    se = numeric(),
    p_value = numeric()
)

# Solo models ----
# Without outcome types (Well-being) ----
for (pred in all_predictors) {
    formula <- paste0("outcome_value ~ ", pred, " + (1 + ", pred, " | outcome_name_id)")
    mdl <- lmer(formula, data = dt)
    s <- summary(mdl)
    coef_table <- s$coefficients
    R2 <- R2Mixed(mdl)
    mixed_results <- rbind(mixed_results, data.table(
        outcome_type = "Well-being",
        model = "solo",
        predictor = pred,
        R2 = R2,
        R2_over = NA_real_,
        beta = if (pred %in% names(fixef(mdl))) fixef(mdl)[pred] else NA_real_,
        se = if (pred %in% rownames(coef_table)) coef_table[pred, "Std. Error"] else NA_real_,
        p_value = if (pred %in% rownames(coef_table)) coef_table[pred, "Pr(>|t|)"] else NA_real_
    ))
}

# With outcome types ----
for (pred in all_predictors) {
    formula <- paste0("outcome_value ~ outcome_type * ", pred, " + (1 + ", pred, " | outcome_name_id)")
    mdl <- lmer(formula, data = dt)
    s <- summary(mdl)
    coef_table <- s$coefficients
    R2 <- R2Mixed(mdl, dt$outcome_type)
    for (otype in outcome_types) {
        coef_name <- paste0("outcome_type", otype, ":", pred)
        mixed_results <- rbind(mixed_results, data.table(
            outcome_type = otype,
            model = "solo",
            predictor = pred,
            R2 = R2[otype],
            R2_over = NA_real_,
            beta = if (coef_name %in% names(fixef(mdl))) fixef(mdl)[coef_name] else fixef(mdl)[pred],
            se = if (coef_name %in% rownames(coef_table)) coef_table[coef_name, "Std. Error"] else coef_table[pred, "Std. Error"],
            p_value = if (coef_name %in% rownames(coef_table)) coef_table[coef_name, "Pr(>|t|)"] else coef_table[pred, "Pr(>|t|)"]
        ))
    }
}

# Over models ----
# Without outcome types (Well-being) ----
# Bench ----
formula <- paste0(
    "outcome_value ~ ", paste(bench, collapse = " + "),
    " + (1 + ", paste(bench, collapse = " + "), " | outcome_name_id)"
)
mdl <- lmer(formula, data = dt)
R2_bench <- as.numeric(R2Mixed(mdl))
for (pred in setdiff(all_predictors, bench)) {
    formula <- paste0(
        "outcome_value ~ ", paste(c(bench, pred), collapse = " + "),
        " + (1 + ", paste(c(bench, pred), collapse = " + "), " | outcome_name_id)"
    )
    mdl <- lmer(formula, data = dt)
    s <- summary(mdl)
    coef_table <- s$coefficients
    R2 <- as.numeric(R2Mixed(mdl))
    mixed_results <- rbind(mixed_results, data.table(
        outcome_type = "Well-being",
        model = "over_bench",
        predictor = pred,
        R2 = R2,
        R2_over = if (!is.na(R2_bench)) R2 - R2_bench else NA_real_,
        beta = if (pred %in% names(fixef(mdl))) fixef(mdl)[pred] else NA_real_,
        se = if (pred %in% rownames(coef_table)) coef_table[pred, "Std. Error"] else NA_real_,
        p_value = if (pred %in% rownames(coef_table)) coef_table[pred, "Pr(>|t|)"] else NA_real_
    ))
}

# Ext ----
formula <- paste0(
    "outcome_value ~ ", paste(ext, collapse = " + "),
    " + (1 + ", paste(ext, collapse = " + "), " | outcome_name_id)"
)
mdl <- lmer(formula, data = dt)
s <- summary(mdl)
coef_table <- s$coefficients
R2_ext <- as.numeric(R2Mixed(mdl))
for (pred in setdiff(all_predictors, ext)) {
    formula <- paste0(
        "outcome_value ~ ", paste(c(ext, pred), collapse = " + "),
        " + (1 + ", paste(c(ext, pred), collapse = " + "), " | outcome_name_id)"
    )
    mdl <- lmer(formula, data = dt)
    s <- summary(mdl)
    coef_table <- s$coefficients
    R2 <- as.numeric(R2Mixed(mdl))
    mixed_results <- rbind(mixed_results, data.table(
        outcome_type = "Well-being",
        model = "over_ext",
        predictor = pred,
        R2 = R2,
        R2_over = if (!is.na(R2_ext)) R2 - R2_ext else NA_real_,
        beta = if (pred %in% names(fixef(mdl))) fixef(mdl)[pred] else NA_real_,
        se = if (pred %in% rownames(coef_table)) coef_table[pred, "Std. Error"] else NA_real_,
        p_value = if (pred %in% rownames(coef_table)) coef_table[pred, "Pr(>|t|)"] else NA_real_
    ))
}

# With outcome_types ----
# Bench ----
formula <- paste0(
    "outcome_value ~ outcome_type * (", paste(bench, collapse = " + "),
    ") + (1 + ", paste(bench, collapse = " + "), " | outcome_name_id)"
)
mdl <- lmer(formula, data = dt)
R2_bench <- R2Mixed(mdl, dt$outcome_type)
for (pred in setdiff(all_predictors, bench)) {
    formula <- paste0(
        "outcome_value ~ outcome_type * (", paste(c(bench, pred), collapse = " + "),
        ") + (1 + ", paste(c(bench, pred), collapse = " + "), " | outcome_name_id)"
    )
    mdl <- lmer(formula, data = dt)
    s <- summary(mdl)
    coef_table <- s$coefficients
    R2 <- R2Mixed(mdl, dt$outcome_type)
    for (otype in outcome_types) {
        coef_name <- paste0("outcome_type", otype, ":", pred)
        mixed_results <- rbind(mixed_results, data.table(
            outcome_type = otype,
            model = "over_bench",
            predictor = pred,
            R2 = R2[otype],
            R2_over = if (!is.null(R2_bench)) R2[otype] - R2_bench[otype] else NA_real_,
            beta = if (coef_name %in% names(fixef(mdl))) fixef(mdl)[coef_name] else fixef(mdl)[pred],
            se = if (coef_name %in% rownames(coef_table)) coef_table[coef_name, "Std. Error"] else coef_table[pred, "Std. Error"],
            p_value = if (coef_name %in% rownames(coef_table)) coef_table[coef_name, "Pr(>|t|)"] else coef_table[pred, "Pr(>|t|)"]
        ))
    }
}

# Ext ----
formula <- paste0(
    "outcome_value ~ outcome_type * (", paste(ext, collapse = " + "),
    ") + (1 + ", paste(ext, collapse = " + "), " | outcome_name_id)"
)
mdl <- lmer(formula, data = dt)
R2_ext <- R2Mixed(mdl, dt$outcome_type)
for (pred in setdiff(all_predictors, ext)) {
    formula <- paste0(
        "outcome_value ~ outcome_type * (", paste(c(ext, pred), collapse = " + "),
        ") + (1 + ", paste(c(ext, pred), collapse = " + "), " | outcome_name_id)"
    )
    mdl <- lmer(formula, data = dt)
    R2 <- R2Mixed(mdl, dt$outcome_type)
    s <- summary(mdl)
    coef_table <- s$coefficients
    for (otype in outcome_types) {
        coef_name <- paste0("outcome_type", otype, ":", pred)
        mixed_results <- rbind(mixed_results, data.table(
            outcome_type = otype,
            model = "over_ext",
            predictor = pred,
            R2 = R2[otype],
            R2_over = if (!is.null(R2_ext)) R2[otype] - R2_ext[otype] else NA_real_,
            beta = if (coef_name %in% names(fixef(mdl))) fixef(mdl)[coef_name] else fixef(mdl)[pred],
            se = if (coef_name %in% rownames(coef_table)) coef_table[coef_name, "Std. Error"] else coef_table[pred, "Std. Error"],
            p_value = if (coef_name %in% rownames(coef_table)) coef_table[coef_name, "Pr(>|t|)"] else coef_table[pred, "Pr(>|t|)"]
        ))
    }
}

# Save results to file
fwrite(mixed_results, "results/R2mixed.csv")
cat("Results saved to: results/R2mixed.csv\n")

# Add R2_plot column
mixed_results[, R2_plot := ifelse(model == "solo", R2, R2_over)]
mixed_results[, R2_plot := as.numeric(R2_plot)]

# Arrange plot outcome type plot order
plot_group_names <- c("Well-being", "Depression", "Borderline", "Satisfaction")

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

# Left join to ensure all valid combinations exist (missing ones become NA)
plot_data <- valid_combinations[mixed_results, on = .(predictor, model, outcome_type)]

# Ensure factors for plotting
plot_data[, model := factor(model, levels = c("solo", "over_bench", "over_ext"))]
plot_data[, predictor := factor(predictor, levels = all_predictors)]
plot_data[, outcome_type := factor(outcome_type, levels = plot_group_names)]

# MATLAB default colors
matlab_cols <- c("solo" = "#0072BD", "over_bench" = "#D95319", "over_ext" = "#EDB120")

# Faceted bar plot for all outcome types with consistent spacing
p <- ggplot(plot_data, aes(x = predictor, y = R2_plot, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.3, preserve = "single"), 
             width = 0.3, na.rm = FALSE) +
    scale_fill_manual(
        values = matlab_cols, name = "Model Type",
        labels = c("Solo", paste0("over ", bench2plot), paste0("over ", ext2plot))
    ) +
    labs(
        title = "R2 - Mixed Models",
        x = "Predictor",
        y = expression(Delta ~ R^2)
    ) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
    ) +
    facet_wrap(~outcome_type)


ggsave("plots/R2mixed.pdf", p, width = 12, height = 7)
cat("Plot saved to: plots/R2mixed.pdf\n")