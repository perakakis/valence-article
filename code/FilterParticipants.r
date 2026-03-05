# FilterParticipants.R: Filter subjects based on compliance and variance thresholds
# Function: FilterParticipants
# Inputs:
#   data: data.table containing PID, day, beep, valence columns
#   Rfull, Cfull: Regression and classification data (optional)
#   minCompliance: Minimum compliance threshold (e.g., 0.75)
#   filterConsecMissing: If TRUE, remove subjects with >= 1 day of
#     consecutive missing valence data (based on modal beeps per day)
# Outputs:
#   List with filtered data and summary stats

FilterParticipants <- function(data, Rfull, minCompliance) {
  library(data.table)

  dt <- as.data.table(data)

  ui <- unique(dt$PID)
  n_subj <- length(ui)

  # Compute modal beeps per day (used for consecutive missing filter)
  beeps_per_day_tab <- dt[, .N, by = .(PID, day)]
  modal_beeps_day <- {
    ux <- unique(beeps_per_day_tab$N)
    ux[which.max(tabulate(match(beeps_per_day_tab$N, ux)))]
  }

  # Compute compliance per subject (vectorised)
  comp_dt <- dt[, .(
    compliance = sum(!is.na(valence)) / .N
  ), by = PID]

  # Flag low compliance
  comp_dt[, ok := compliance > minCompliance]
  deleted_low <- comp_dt[ok == FALSE, .N]

  # # Consecutive missing filter
  # deleted_consec <- 0L
  # if (filterConsecMissing) {
  #   consec_dt <- dt[, {
  #     is_na_vec <- is.na(valence)
  #     if (any(is_na_vec)) {
  #       na_runs <- rle(is_na_vec)
  #       max_consec <- max(na_runs$lengths[na_runs$values])
  #     } else {
  #       max_consec <- 0L
  #     }
  #     .(max_consec_na = max_consec)
  #   }, by = PID]

  #   comp_dt <- merge(comp_dt, consec_dt, by = "PID")
  #   deleted_consec <- comp_dt[ok == TRUE & max_consec_na >= modal_beeps_day, .N]
  #   comp_dt[ok == TRUE & max_consec_na >= modal_beeps_day, ok := FALSE]
  # }

  # Apply filter
  keep_pids <- comp_dt[ok == TRUE, PID]

  if (length(keep_pids) == 0) {
    stop("No subjects remain after filtering")
  }

  dt_new <- dt[PID %in% keep_pids]

  # Filter regression outcomes (intersect: keep only subjects present in both)
  deleted_no_reg <- 0L
  if (!is.null(Rfull)) {
    Rfull_dt <- as.data.table(Rfull)
    setnames(Rfull_dt, 1, "PID")
    RfullNew <- Rfull_dt[PID %in% keep_pids]
    missing_reg <- setdiff(keep_pids, RfullNew$PID)
    deleted_no_reg <- length(missing_reg)
    if (deleted_no_reg > 0L) {
      message(sprintf("  %d subject(s) removed: present in data but missing from reg.csv",
                      deleted_no_reg))
      keep_pids <- intersect(keep_pids, RfullNew$PID)
      dt_new <- dt_new[PID %in% keep_pids]
      RfullNew <- RfullNew[PID %in% keep_pids]
    }
  } else {
    RfullNew <- NULL
  }

  # # Filter classification outcomes (intersect: keep only subjects present in both)
  # deleted_no_class <- 0L
  # if (!is.null(Cfull)) {
  #   Cfull_dt <- as.data.table(Cfull)
  #   setnames(Cfull_dt, 1, "PID")
  #   CfullNew <- Cfull_dt[PID %in% keep_pids]
  #   missing_class <- setdiff(keep_pids, CfullNew$PID)
  #   deleted_no_class <- length(missing_class)
  #   if (deleted_no_class > 0L) {
  #     message(sprintf("  %d subject(s) removed: present in data but missing from class.csv",
  #                     deleted_no_class))
  #     keep_pids <- intersect(keep_pids, CfullNew$PID)
  #     dt_new <- dt_new[PID %in% keep_pids]
  #     RfullNew <- if (!is.null(RfullNew)) RfullNew[PID %in% keep_pids] else NULL
  #     CfullNew <- CfullNew[PID %in% keep_pids]
  #   }
  # } else {
  #   CfullNew <- NULL
  # }

  # if (length(keep_pids) == 0) {
  #   stop("No subjects remain after filtering (including outcome matching)")
  # }

  # Summary stats
  ok_compliance <- comp_dt[ok == TRUE, compliance]
  dataPointsInfo <- c(
    meanCompliance = mean(ok_compliance, na.rm = TRUE),
    meanBeepM = mean(dt_new[, .N, by = PID]$N)
  )

  deletedInfo <- c(
    totalSubj = n_subj,
    subjKept = length(keep_pids),
    subjDeleted = n_subj - length(keep_pids),
    lowCompliance = deleted_low,
    # consecMissing = deleted_consec,
    noRegOutcome = deleted_no_reg
    # noClassOutcome = deleted_no_class
  )

  # Median days per subject (more representative than global unique days)
  median_days <- median(dt_new[, .(n_days = uniqueN(day)), by = PID]$n_days)

  return(list(
    dataNew = dt_new,
    RfullNew = RfullNew,
    # CfullNew = CfullNew,
    deletedInfo = deletedInfo,
    dataPointsInfo = dataPointsInfo,
    medianDays = median_days,
    modalBeepsPerDay = modal_beeps_day,
    keptPIDs = keep_pids
  ))
}
