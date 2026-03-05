# Generic transformation script for EMA and survey data
# Transform both EMA and survey files together to ensure participant alignment

# ============================================================================
# CONFIGURATION - Edit this section to switch between datasets
# ============================================================================

DATASET <- "postcovid1"  # Options: "german", "postcovid1", "postcovid2", etc.

# Dataset-specific configurations
if (DATASET == "german") {
  # German dataset
  ema_input <- "data/00 german_EMA.csv"
  survey_input <- "data/00 german_survey.csv"
  ema_output <- "data/00 german.csv"
  reg_output <- "data/00 german reg.csv"

  # Column mappings for EMA
  valence_col <- "Core-Affect-Valence"
  id_col <- "VID"
  day_col <- "Day"
  session_col <- "Session"

  # Valence transformation
  valence_offset <- -50  # Shift from 0-100 to -50-50
  valence_split_point <- 50  # Split point for PA/NA

  # Survey outcome computation function
  compute_outcomes <- function(df) {
    df$DASSd <- df$DASS.Depression
    df$DASSa <- df$DASS.Anxiety
    df$BRS <- df$BriefResilience
    df$AAQ <- df$AAQ
    df$SWL <- df$SatisfactLife
    df <- df[c("PID", "DASSd", "DASSa", "AAQ", "BRS", "SWL")]
    return(df)
  }

} else if (DATASET == "postcovid1") {
  # Postcovid1 dataset
  dataset_prefix <- "01 postcovid1"
  ema_input <- "data/01 postcovid1_ema.csv"
  survey_input <- "data/01 postcovid1_survey.csv"
  ema_output <- "data/01 postcovid1.csv"
  reg_output <- "data/01 postcovid1 reg.csv"

  # Column mappings for EMA
  valence_col <- "valence"
  id_col <- "participant"
  day_col <- "study_day"
  session_col <- "occasion"

  # Valence transformation
  valence_offset <- 50  # Shift from -50-50 to 0-100
  valence_split_point <- 0  # Split point for PA/NA (before offset)

  # Survey outcome computation function
  compute_outcomes <- function(df) {
    df$PHQ <- rowSums(df[, c("PHQ_1", "PHQ_2", "PHQ_3", "PHQ_4", "PHQ_5", "PHQ_6", "PHQ_7", "PHQ_8", "PHQ_9")])
    df$GAD <- rowSums(df[, c("GAD_1", "GAD_2", "GAD_3", "GAD_4", "GAD_5", "GAD_6", "GAD_7")])
    df$BRS <- rowMeans(cbind(df$BRS_1, df$BRS_3, df$BRS_5, 6 - df$BRS_2, 6 - df$BRS_4, 6 - df$BRS_6))
    df$AAQ <- rowSums(df[, paste0("AAQ_", 1:7)])
    df$FS <- rowSums(df[, paste0("FS_", 1:8)])
    df$SWL <- df$GLS
    df <- df[c("PID", "PHQ", "GAD", "BRS", "AAQ", "FS", "SWL")]
    return(df)
  }

} else {
  stop("Unknown dataset: ", DATASET)
}

# ============================================================================
# MAIN TRANSFORMATION - Do not edit below this line
# ============================================================================

# Load required libraries
library(readr)
library(dplyr)

cat("==============================================\n")
cat("Transforming dataset:", DATASET, "\n")
cat("==============================================\n\n")

# ============================================================================
# PART 1: Transform EMA data
# ============================================================================

cat("PART 1: Transforming EMA data\n")
cat("----------------------------------------------\n")

# Read EMA data
cat("Reading EMA data from:", ema_input, "\n")
ema_data <- read_csv(ema_input, show_col_types = FALSE)

# Rename columns based on mappings
names(ema_data)[names(ema_data) == valence_col] <- "valence"
names(ema_data)[names(ema_data) == id_col] <- "PID"
names(ema_data)[names(ema_data) == day_col] <- "UNIT"
names(ema_data)[names(ema_data) == session_col] <- "OCCASION"

cat("Original EMA dimensions:", nrow(ema_data), "rows,", ncol(ema_data), "columns\n")

# Keep only relevant columns
ema_data <- ema_data[c("PID", "UNIT", "OCCASION", "valence")]

# Get unique participants (this establishes the canonical order)
unique_participants <- unique(ema_data$PID)
cat("Number of unique participants:", length(unique_participants), "\n")

# Keep original PID (no transformation to numeric)
# unique_participants maintains consistent ordering across EMA and survey

# Create PAf (Positive Affect) and NAf (Negative Affect)
# These are computed BEFORE applying the offset transformation
# Note: Values > split point go to PAf, values < split point go to NAf
# Values equal to split point (neutral) are excluded from both (set to NA)
ema_data$PAf <- ifelse(!is.na(ema_data$valence) & ema_data$valence > valence_split_point,
                       abs(ema_data$valence - valence_split_point), NA)

ema_data$NAf <- ifelse(!is.na(ema_data$valence) & ema_data$valence < valence_split_point,
                       abs(ema_data$valence - valence_split_point), NA)

# Apply valence offset transformation (converts to 0-100 scale)
ema_data$valence <- ema_data$valence + valence_offset

# Sort by PID, UNIT, OCCASION
ema_data <- ema_data[order(ema_data$PID, ema_data$UNIT, ema_data$OCCASION), ]

cat("Final EMA dimensions:", nrow(ema_data), "rows,", ncol(ema_data), "columns\n")

# Rename columns to match INITIALISE.R expectations right before export
ema_export <- ema_data
names(ema_export)[names(ema_export) == "valence"] <- "Happy"
names(ema_export)[names(ema_export) == "PAf"] <- "PA"
names(ema_export)[names(ema_export) == "NAf"] <- "NA"

# Export EMA data
cat("Writing EMA data to:", ema_output, "\n")
write.table(ema_export, file = ema_output, sep = ";", row.names = FALSE, quote = FALSE, na = "")
cat("EMA transformation complete!\n\n")

# ============================================================================
# PART 2: Transform survey data
# ============================================================================

cat("PART 2: Transforming survey data\n")
cat("----------------------------------------------\n")

# Read survey data
cat("Reading survey data from:", survey_input, "\n")
survey_data <- read_csv(survey_input, show_col_types = FALSE)

cat("Original survey dimensions:", nrow(survey_data), "rows,", ncol(survey_data), "columns\n")

# Rename participant column to match
if (id_col %in% names(survey_data)) {
  names(survey_data)[names(survey_data) == id_col] <- "participant_original"
} else {
  names(survey_data)[names(survey_data) == "participant"] <- "participant_original"
}

# Use original PID (match ordering from EMA via unique_participants)
survey_data$PID <- factor(survey_data$participant_original, levels = unique_participants)
survey_data$PID <- as.character(survey_data$PID)  # Convert back to character

# Compute outcomes using dataset-specific function
survey_data <- compute_outcomes(survey_data)

cat("Final survey dimensions:", nrow(survey_data), "rows,", ncol(survey_data), "columns\n")

# Export survey data
cat("Writing survey data to:", reg_output, "\n")
write.csv(survey_data, file = reg_output, row.names = FALSE)
cat("Survey transformation complete!\n\n")

# ============================================================================
# VERIFICATION & QUALITY CHECKS
# ============================================================================

cat("==============================================\n")
cat("VERIFICATION & QUALITY CHECKS\n")
cat("==============================================\n\n")

# Check 1: Participant alignment
cat("CHECK 1: Participant Alignment\n")
cat("----------------------------------------------\n")
ema_pids <- sort(unique(ema_data$PID))
survey_pids <- sort(unique(survey_data$PID))
cat("EMA participants:", length(ema_pids), "unique PIDs\n")
cat("Survey participants:", length(survey_pids), "unique PIDs\n")

if (identical(ema_pids, survey_pids)) {
  cat("✓ PASS: Participant IDs match perfectly between EMA and survey\n\n")
} else {
  cat("✗ FAIL: Participant IDs DO NOT match\n")
  cat("  Missing in survey:", setdiff(ema_pids, survey_pids), "\n")
  cat("  Missing in EMA:", setdiff(survey_pids, ema_pids), "\n\n")
}

# Check 2: Valence range verification
cat("CHECK 2: Valence Transformation\n")
cat("----------------------------------------------\n")
valence_stats <- summary(ema_data$valence[!is.na(ema_data$valence)])
cat("Transformed valence range: [", valence_stats[1], ",", valence_stats[6], "]\n", sep = "")
cat("Expected range: [0, 100]\n")
if (valence_stats[1] >= 0 && valence_stats[6] <= 100) {
  cat("✓ PASS: Valence is in expected range [0, 100]\n\n")
} else {
  cat("✗ FAIL: Valence is outside expected range\n\n")
}

# Check 3: PAf and NAf verification
cat("CHECK 3: Positive/Negative Affect (PAf/NAf)\n")
cat("----------------------------------------------\n")
paf_count <- sum(!is.na(ema_data$PAf))
naf_count <- sum(!is.na(ema_data$NAf))
both_count <- sum(!is.na(ema_data$PAf) & !is.na(ema_data$NAf))
total_valid <- sum(!is.na(ema_data$valence))
cat("Total valid valence entries:", total_valid, "\n")
cat("PAf (Positive Affect) entries:", paf_count, sprintf("(%.1f%%)\n", 100*paf_count/total_valid))
cat("NAf (Negative Affect) entries:", naf_count, sprintf("(%.1f%%)\n", 100*naf_count/total_valid))
cat("Both PAf and NAf (should be 0):", both_count, "\n")
if (both_count == 0 && (paf_count + naf_count) == total_valid) {
  cat("✓ PASS: PAf/NAf are mutually exclusive and complete\n")
} else {
  cat("✗ WARNING: PAf/NAf coverage issue detected\n")
}
cat("PAf range: [", min(ema_data$PAf, na.rm = TRUE), ",", max(ema_data$PAf, na.rm = TRUE), "]\n", sep = "")
cat("NAf range: [", min(ema_data$NAf, na.rm = TRUE), ",", max(ema_data$NAf, na.rm = TRUE), "]\n\n", sep = "")

# Check 4: Survey outcome completeness
cat("CHECK 4: Survey Outcome Variables\n")
cat("----------------------------------------------\n")
outcome_cols <- setdiff(names(survey_data), "PID")
for (col in outcome_cols) {
  missing <- sum(is.na(survey_data[[col]]))
  missing_pct <- 100 * missing / nrow(survey_data)
  cat(sprintf("%-10s: %3d missing (%.1f%%) | Range: [%.2f, %.2f] | Mean: %.2f\n",
              col, missing, missing_pct,
              min(survey_data[[col]], na.rm = TRUE),
              max(survey_data[[col]], na.rm = TRUE),
              mean(survey_data[[col]], na.rm = TRUE)))
}
cat("\n")

# Check 5: Sample data verification
cat("CHECK 5: Sample Data Verification\n")
cat("----------------------------------------------\n")
cat("First 3 EMA rows for PID=1:\n")
sample_ema <- head(ema_data[ema_data$PID == 1, ], 3)
print(sample_ema)
cat("\nFirst survey row for PID=1:\n")
sample_survey <- head(survey_data[survey_data$PID == 1, ], 1)
print(sample_survey)
cat("\n")

# Check 6: Data completeness per participant
cat("CHECK 6: EMA Data Completeness per Participant\n")
cat("----------------------------------------------\n")
ema_counts <- aggregate(valence ~ PID, data = ema_data, FUN = function(x) sum(!is.na(x)))
names(ema_counts)[2] <- "valid_entries"
cat("Average valid entries per participant:", round(mean(ema_counts$valid_entries), 1), "\n")
cat("Min valid entries:", min(ema_counts$valid_entries), "\n")
cat("Max valid entries:", max(ema_counts$valid_entries), "\n")
low_compliance <- ema_counts$PID[ema_counts$valid_entries < 10]
if (length(low_compliance) > 0) {
  cat("⚠ WARNING:", length(low_compliance), "participants with < 10 valid entries\n")
  cat("  PIDs:", paste(head(low_compliance, 10), collapse = ", "),
      ifelse(length(low_compliance) > 10, "...", ""), "\n")
} else {
  cat("✓ All participants have >= 10 valid entries\n")
}

cat("\n==============================================\n")
cat("GENERATING VERIFICATION PLOTS\n")
cat("==============================================\n\n")

# Create plots directory if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Read original EMA data for plotting
original_ema <- read.csv(ema_input)
original_ema$PID <- factor(original_ema[[id_col]], levels = unique_participants)
original_ema$PID <- as.character(original_ema$PID)

# Plot 1: Valence distribution comparison (original vs transformed)
pdf(paste0("plots/", dataset_prefix, " valence_transformation.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

# Original valence
hist(original_ema[[valence_col]],
     main = paste(DATASET, "- Original Valence Distribution"),
     xlab = "Valence (original scale)",
     col = "lightblue",
     breaks = 30)
abline(v = valence_split_point, col = "red", lwd = 2, lty = 2)
text(valence_split_point, par("usr")[4] * 0.9,
     paste("Split point =", valence_split_point),
     pos = 4, col = "red")

# Transformed valence
hist(ema_data$valence,
     main = paste(DATASET, "- Transformed Valence Distribution"),
     xlab = "Valence (0-100 scale)",
     col = "lightgreen",
     breaks = 30)
abline(v = 50, col = "red", lwd = 2, lty = 2)
text(50, par("usr")[4] * 0.9, "Split point = 50", pos = 4, col = "red")

dev.off()
cat("✓ Saved: plots/", dataset_prefix, " valence_transformation.pdf\n", sep = "")

# Plot 2: Valence trajectories per subject - Original vs Transformed (side by side)
pdf(paste0("plots/", dataset_prefix, " valence_by_subject.pdf"), width = 16, height = 10)
par(mfrow = c(3, 6), mar = c(3, 3, 2.5, 1), oma = c(0, 0, 2, 0))

# Select 9 subjects with most complete data
subject_completeness <- aggregate(valence ~ PID, data = ema_data,
                                  FUN = function(x) sum(!is.na(x)))
top_subjects <- head(subject_completeness[order(-subject_completeness$valence), "PID"], 9)

for (pid in top_subjects) {
  # Original data
  orig_subj <- original_ema[original_ema$PID == pid, ]
  orig_vals <- orig_subj[[valence_col]]

  # Transformed data
  trans_subj <- ema_data[ema_data$PID == pid, ]
  trans_vals <- trans_subj$valence

  # Create time index
  time_idx <- seq_along(orig_vals)

  # Plot original valence
  plot(time_idx, orig_vals,
       type = "l", col = "blue", lwd = 2,
       main = paste("Subject", pid, "- Original"),
       xlab = "Occasion",
       ylab = "Valence",
       ylim = range(orig_vals, na.rm = TRUE) + c(-5, 5))
  abline(h = valence_split_point, col = "red", lty = 2)

  # Plot transformed valence (actual 0-100 scale)
  plot(time_idx, trans_vals,
       type = "l", col = "darkgreen", lwd = 2,
       main = paste("Subject", pid, "- Transformed"),
       xlab = "Occasion",
       ylab = "Valence",
       ylim = c(0, 100))
  abline(h = 50, col = "red", lty = 2)
}

mtext("Valence Time Series: Original vs Transformed Scales", outer = TRUE, cex = 1.2, font = 2)

dev.off()
cat("✓ Saved: plots/", dataset_prefix, " valence_by_subject.pdf\n", sep = "")

# Plot 3: PAf vs NAf distribution
pdf(paste0("plots/", dataset_prefix, " PAf_NAf_distribution.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))

# PAf distribution
hist(ema_data$PAf,
     main = paste(DATASET, "- Positive Affect (PAf) Distribution"),
     xlab = "PAf (distance from split point)",
     col = "lightcoral",
     breaks = 30)

# NAf distribution
hist(ema_data$NAf,
     main = paste(DATASET, "- Negative Affect (NAf) Distribution"),
     xlab = "NAf (distance from split point)",
     col = "lightyellow",
     breaks = 30)

dev.off()
cat("✓ Saved: plots/", dataset_prefix, " PAf_NAf_distribution.pdf\n", sep = "")

# Plot 4: Subject-level summaries
pdf(paste0("plots/", dataset_prefix, " subject_summaries.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2))

# Mean valence per subject (original)
orig_means <- aggregate(original_ema[[valence_col]] ~ PID,
                        data = original_ema, FUN = mean, na.rm = TRUE)
barplot(orig_means[, 2],
        names.arg = orig_means$PID,
        main = "Mean Valence per Subject (Original Scale)",
        xlab = "Subject ID",
        ylab = "Mean Valence",
        col = "skyblue",
        las = 2,
        cex.names = 0.6)
abline(h = valence_split_point, col = "red", lty = 2)

# Mean valence per subject (transformed)
trans_means <- aggregate(valence ~ PID, data = ema_data, FUN = mean, na.rm = TRUE)
barplot(trans_means$valence,
        names.arg = trans_means$PID,
        main = "Mean Valence per Subject (Transformed 0-100)",
        xlab = "Subject ID",
        ylab = "Mean Valence",
        col = "lightgreen",
        las = 2,
        cex.names = 0.6)
abline(h = 50, col = "red", lty = 2)

# PAf percentage per subject
paf_counts <- aggregate(cbind(PAf = PAf, NAf = NAf) ~ PID, data = ema_data,
                        FUN = function(x) sum(!is.na(x)), na.action = na.pass)
paf_counts$pct <- 100 * paf_counts$PAf / (paf_counts$PAf + paf_counts$NAf)
barplot(paf_counts$pct,
        names.arg = paf_counts$PID,
        main = "% Positive Affect per Subject",
        xlab = "Subject ID",
        ylab = "% Measurements with PAf",
        col = "coral",
        las = 2,
        cex.names = 0.6)
abline(h = 50, col = "red", lty = 2)

# Valid measurements per subject
valid_counts <- aggregate(valence ~ PID, data = ema_data,
                          FUN = function(x) sum(!is.na(x)))
barplot(valid_counts$valence,
        names.arg = valid_counts$PID,
        main = "Valid Measurements per Subject",
        xlab = "Subject ID",
        ylab = "Count",
        col = "lightblue",
        las = 2,
        cex.names = 0.6)

dev.off()
cat("✓ Saved: plots/", dataset_prefix, " subject_summaries.pdf\n", sep = "")

# Plot 5: Survey outcome verification - Original vs Transformed participant matching
cat("\nGenerating survey outcome verification plots...\n")

# Read original survey data
original_survey <- read_csv(survey_input, show_col_types = FALSE)

# Add participant_original column
if (id_col %in% names(original_survey)) {
  original_survey$participant_original <- original_survey[[id_col]]
} else {
  original_survey$participant_original <- original_survey$participant
}

original_survey$PID <- factor(original_survey$participant_original, levels = unique_participants)
original_survey$PID <- as.character(original_survey$PID)

# Get outcome columns from transformed data
outcome_cols <- setdiff(names(survey_data), "PID")
n_outcomes <- length(outcome_cols)

# Calculate number of rows needed for subplot grid
n_rows <- ceiling(sqrt(n_outcomes))
n_cols <- ceiling(n_outcomes / n_rows)

pdf(paste0("plots/", dataset_prefix, " survey_verification.pdf"), width = 4 * n_cols, height = 4 * n_rows)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 2))

# For each outcome, plot original vs transformed values
for (outcome in outcome_cols) {
  # Get values from transformed survey
  trans_values <- survey_data[[outcome]]

  # Compute original values using the same function
  original_computed <- compute_outcomes(original_survey)
  orig_values <- original_computed[[outcome]]

  # Create scatter plot
  plot(orig_values, trans_values,
       xlab = "Original Survey Data",
       ylab = "Transformed Survey Data",
       main = paste(outcome, "- Participant Matching"),
       pch = 19,
       col = rgb(0, 0, 1, 0.6))

  # Add 1:1 reference line
  abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

  # Add correlation
  cor_val <- cor(orig_values, trans_values, use = "complete.obs")
  text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
       y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
       labels = sprintf("r = %.4f", cor_val),
       adj = c(0, 1),
       cex = 0.9,
       col = ifelse(abs(cor_val - 1) < 0.0001, "darkgreen", "red"))

  # Check for perfect match
  if (all(abs(orig_values - trans_values) < 1e-10, na.rm = TRUE)) {
    text(x = mean(par("usr")[1:2]),
         y = par("usr")[3] + 0.1 * diff(par("usr")[3:4]),
         labels = "✓ PERFECT MATCH",
         col = "darkgreen",
         font = 2,
         cex = 0.8)
  }
}

dev.off()
cat("✓ Saved: plots/", dataset_prefix, " survey_verification.pdf\n", sep = "")

# Additional text verification
cat("\nCHECK 7: Survey Data Participant Matching\n")
cat("----------------------------------------------\n")
for (outcome in outcome_cols) {
  trans_values <- survey_data[[outcome]]
  original_computed <- compute_outcomes(original_survey)
  orig_values <- original_computed[[outcome]]

  max_diff <- max(abs(orig_values - trans_values), na.rm = TRUE)
  cor_val <- cor(orig_values, trans_values, use = "complete.obs")

  if (max_diff < 1e-10 && abs(cor_val - 1) < 1e-10) {
    cat(sprintf("%-10s: ✓ PERFECT MATCH (max diff = %.2e, r = %.6f)\n", outcome, max_diff, cor_val))
  } else {
    cat(sprintf("%-10s: ✗ MISMATCH (max diff = %.4f, r = %.6f)\n", outcome, max_diff, cor_val))
  }
}
cat("\n")

cat("\n==============================================\n")
cat("TRANSFORMATION COMPLETE FOR:", DATASET, "\n")
cat("==============================================\n")
