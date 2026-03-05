# HealthCheck.R: Diagnostic checks on ESM dataset structure
# Checks: day count consistency, beep count consistency,
#          OCCASION sequentiality, BEEP sequentiality
library(data.table)

# Auto-discover main data CSVs (exclude reg/cov/val/class/bound/cat files)
all_csvs <- list.files("data", pattern = "\\.csv$", full.names = TRUE)
suffix_pattern <- "(reg|cov|val|class|bound|cat)\\.csv$"
main_files <- all_csvs[!grepl(suffix_pattern, all_csvs)]

cat(sprintf("Found %d main data files\n\n", length(main_files)))

# Helper: statistical mode (most frequent value)
stat_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Store summary across datasets
summary_rows <- list()
# Collect abnormal findings for final report
abnormal_findings <- list()

for (f in main_files) {
  dataset_name <- gsub("^data/", "", gsub("\\.csv$", "", f))

  dt <- fread(f, sep = "auto", header = TRUE, na.strings = c(""),
              colClasses = list(character = "PID"), check.names = FALSE)

  # Check required columns exist
  required <- c("PID", "UNIT", "BEEP", "OCCASION")
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0) {
    next
  }

  n_pids <- length(unique(dt$PID))

  # ------------------------------------------------------------------
  # 1. Day count per PID
  # ------------------------------------------------------------------
  days_per_pid <- dt[, .(n_days = length(unique(UNIT))), by = PID]
  modal_days <- stat_mode(days_per_pid$n_days)
  day_min <- min(days_per_pid$n_days)
  day_max <- max(days_per_pid$n_days)
  day_sd <- sd(days_per_pid$n_days)
  day_mean <- mean(days_per_pid$n_days)
  pids_off_days <- sum(days_per_pid$n_days != modal_days)
  pct_off_days <- round(100 * pids_off_days / n_pids, 1)

  # ------------------------------------------------------------------
  # 2. Beep count per UNIT per PID
  # ------------------------------------------------------------------
  beeps_per_unit <- dt[, .(n_beeps = .N), by = .(PID, UNIT)]
  modal_beeps <- stat_mode(beeps_per_unit$n_beeps)
  beep_min <- min(beeps_per_unit$n_beeps)
  beep_max <- max(beeps_per_unit$n_beeps)
  beep_sd <- sd(beeps_per_unit$n_beeps)
  beep_mean <- mean(beeps_per_unit$n_beeps)
  units_off_beeps <- sum(beeps_per_unit$n_beeps != modal_beeps)

  # ------------------------------------------------------------------
  # 3. OCCASION sequentiality
  # ------------------------------------------------------------------
  occasion_check <- dt[, .(
    is_sequential = identical(as.integer(OCCASION), seq_len(.N))
  ), by = PID]
  pids_bad_occasion <- sum(!occasion_check$is_sequential)

  # ------------------------------------------------------------------
  # 4. BEEP sequentiality
  # ------------------------------------------------------------------
  beep_seq_check <- dt[, .(
    is_sequential = identical(as.integer(BEEP), seq_len(.N))
  ), by = .(PID, UNIT)]
  units_bad_beep_seq <- sum(!beep_seq_check$is_sequential)

  # ------------------------------------------------------------------
  # 5. Detect non-contiguous UNIT blocks (same UNIT appears in separate runs)
  # ------------------------------------------------------------------
  unit_runs <- dt[, .(unit_changed = UNIT != shift(UNIT, 1, type = "lag", fill = UNIT[1])), by = PID]
  # Count unique UNITs vs number of UNIT runs per PID
  unit_run_count <- dt[, {
    rle_units <- rle(as.character(UNIT))
    .(n_unique_units = length(unique(UNIT)), n_runs = length(rle_units$values))
  }, by = PID]
  pids_noncontig <- sum(unit_run_count$n_runs > unit_run_count$n_unique_units)

  # ------------------------------------------------------------------
  # 6. Check for UNIT = 0
  # ------------------------------------------------------------------
  has_unit_zero <- any(dt$UNIT == 0, na.rm = TRUE)

  # Store summary row
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    dataset = dataset_name,
    n_pids = n_pids,
    n_rows = nrow(dt),
    modal_days = modal_days,
    day_mean = round(day_mean, 2),
    day_sd = round(day_sd, 2),
    day_min = day_min,
    day_max = day_max,
    pids_off_days = pids_off_days,
    pct_off_days = pct_off_days,
    modal_beeps = modal_beeps,
    beep_mean = round(beep_mean, 2),
    beep_sd = round(beep_sd, 2),
    beep_min = beep_min,
    beep_max = beep_max,
    units_off_beeps = units_off_beeps,
    pids_bad_occasion = pids_bad_occasion,
    units_bad_beep_seq = units_bad_beep_seq,
    pids_noncontig_units = pids_noncontig,
    has_unit_zero = has_unit_zero
  )

  # ------------------------------------------------------------------
  # Detect abnormal patterns (beyond normal ESM missingness)
  # Normal: +-1 day from mode, fewer beeps in some days (missed notifications)
  # Abnormal: large day range, UNIT=0, non-contiguous UNITs, beeps > mode
  # ------------------------------------------------------------------
  findings <- character()

  # Day range too wide (more than +-1 from mode)
  if (day_max - day_min > 2 || day_max > modal_days + 1) {
    findings <- c(findings, sprintf(
      "Day count range %d-%d (mode %d): %d PIDs (%.0f%%) deviate",
      day_min, day_max, modal_days, pids_off_days, pct_off_days))
  }

  # UNIT = 0
  if (has_unit_zero) {
    findings <- c(findings, "UNIT = 0 present (day indexing starts at 0)")
  }

  # Non-contiguous UNITs (same UNIT appears in separate blocks within a PID)
  if (pids_noncontig > 0) {
    findings <- c(findings, sprintf(
      "%d PIDs have non-contiguous UNITs (same day appears in separate blocks)",
      pids_noncontig))
  }

  # Beeps above modal (extra beeps, not just missed ones)
  if (beep_max > modal_beeps) {
    above_modal <- beeps_per_unit[n_beeps > modal_beeps]
    findings <- c(findings, sprintf(
      "%d PID-day combos have MORE beeps than mode (%d), max = %d",
      nrow(above_modal), modal_beeps, beep_max))
  }

  # Non-sequential OCCASION
  if (pids_bad_occasion > 0) {
    findings <- c(findings, sprintf(
      "%d PIDs have non-sequential OCCASION", pids_bad_occasion))
  }

  # Non-sequential BEEP (in pre-existing BEEP columns)
  if (units_bad_beep_seq > 0) {
    findings <- c(findings, sprintf(
      "%d PID-day combos have non-sequential BEEP values", units_bad_beep_seq))
  }

  if (length(findings) > 0) {
    abnormal_findings[[dataset_name]] <- findings
  }
}

# ------------------------------------------------------------------
# Save CSV
# ------------------------------------------------------------------
summary_dt <- rbindlist(summary_rows)
write.csv(summary_dt, "results/healthcheck.csv", row.names = FALSE)
cat("Summary saved to: results/healthcheck.csv\n\n")

# ------------------------------------------------------------------
# Console: only abnormal patterns
# ------------------------------------------------------------------
if (length(abnormal_findings) == 0) {
  cat("All datasets look clean.\n")
} else {
  cat(sprintf("Abnormal patterns found in %d/%d datasets:\n\n",
              length(abnormal_findings), nrow(summary_dt)))
  for (ds in names(abnormal_findings)) {
    cat(sprintf("  %s:\n", ds))
    for (finding in abnormal_findings[[ds]]) {
      cat(sprintf("    - %s\n", finding))
    }
    cat("\n")
  }
}
