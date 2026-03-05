# MEAN_SD.R: Compute mean and standard deviation per subject for PA, NA
# Function: MEAN_SD
# Inputs:
#   data: Data.table with columns subj, beep, PosA, NegA, Valence
#   stat_vars: Vector of column names to compute stats (e.g., c("PA", "NA"))
# Outputs:
#   Data.table with columns subj, M_<var>, SD_<var>, L_<var> (e.g., M_PA, SD_PA, L_PA)

MEAN_SD <- function(data, stat_vars) {
  library(data.table)

  # Validate inputs
  if (!all(stat_vars %in% names(data))) {
    stop("Missing stat_vars in data: ", paste(setdiff(stat_vars, names(data)), collapse = ", "))
  }
  if (!"PID" %in% names(data)) {
    stop("PID column missing in data")
  }

  # Initialize result
  result <- data.table(PID = unique(data$PID))

  # Compute mean, SD, and counts for each variable
  for (var in stat_vars) {
    # Group by PID, compute mean, SD, and count of non-NA
    stats <- data[, .(
      mean_val = mean(get(var), na.rm = TRUE),
      sd_val = sd(get(var), na.rm = TRUE),
      count = sum(!is.na(get(var)))
    ),
    by = PID
    ]

    # Set NA for insufficient data
    stats[count == 0, mean_val := NA_real_]
    stats[count < 2, sd_val := NA_real_]

    # Rename columns (omit count from result)
    setnames(
      stats, c("mean_val", "sd_val"),
      c(paste0("M_", var), paste0("SD_", var))
    )

    # Merge with result (omit count)
    result <- merge(result, stats[, .(
      PID, get(paste0("M_", var)),
      get(paste0("SD_", var))
    )],
    by = "PID"
    )
    setnames(
      result, c(paste0("V", 2:3)),
      c(paste0("M_", var), paste0("SD_", var))
    )
  }

  # Check for subjects with all NA
  if (all(is.na(result[, paste0("M_", stat_vars), with = FALSE])) &&
    all(is.na(result[, paste0("SD_", stat_vars), with = FALSE]))) {
    warning("All subjects have NA means and SDs for all variables")
  }

  return(result)
}