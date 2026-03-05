# ASR.R: Compute P2N and N2P transitions per subject for Valence
# Function: ASR
# Inputs:
#   data: Data.table with columns PID, day (optional), beep, and the variable specified in asr_var
#   asr_var: Column name for transitions (e.g., "valence")
#   center_value: Value to subtract for zero-centering (e.g., 50)
# Outputs:
#   Data.table with columns PID, P2N, N2P
#
# Note: Data should be pre-sorted by day and beep within each PID, or include a 'day' column for proper sorting

ASR <- function(data, asr_var, center_value = 50) {
  library(data.table)

  # Validate inputs
  if (!asr_var %in% names(data)) {
    stop("Missing asr_var in data: ", asr_var)
  }
  if (!"PID" %in% names(data)) {
    stop("PID column missing in data")
  }
  if (!"beep" %in% names(data)) {
    stop("beep column missing in data")
  }

  # Copy data and center Valence
  dt <- copy(data)
  dt[, (asr_var) := get(asr_var) - center_value]

  # Initialize result
  result <- data.table(PID = unique(dt$PID))
  result[, c("P2N", "N2P") := .(NA_real_, NA_real_)]

  # Compute transitions per subject
  for (s in result$PID) {
    # Get subject data (assumes data is already sorted by PID, day, beep)
    PID_dt <- dt[PID == s]
    val <- PID_dt[[asr_var]]

    # Remove NAs
    # val <- na.omit(val)

    # Count non-NA observations
    n_valid <- sum(!is.na(val))

    # Skip if insufficient data
    if (n_valid < 2) next

    # Count positive and negative observations (exclude y = 0)
    # Keep NAs in the sequence to match Spanish method
    n_pos <- sum(val > 0, na.rm = TRUE)
    n_neg <- sum(val < 0, na.rm = TRUE)

    # Compute transitions using diff(sign(y)) - Spanish method
    # This keeps NAs in sequence, so only consecutive observations are compared
    p2n_indices <- which(diff(sign(val)) == -2)
    n2p_indices <- which(diff(sign(val)) == 2)

    p2n_trans <- length(p2n_indices)
    n2p_trans <- length(n2p_indices)

    # Normalize transitions - Spanish method
    # Only assign value if transitions exist, otherwise leave as NA
    if (p2n_trans > 0) {
      p2n <- p2n_trans / n_pos
    } else {
      p2n <- NA_real_  # No transitions or no positive observations - leave as NA
    }

    if (n2p_trans > 0) {
      n2p <- n2p_trans / n_neg
    } else {
      n2p <- NA_real_  # No transitions or no negative observations - leave as NA
    }

    result[PID == s, c("P2N", "N2P") := .(p2n, n2p)]
  }

  # Check for subjects with all NA
  if (all(is.na(result[, .(P2N, N2P)]))) {
    warning("All subjects have NA P2N and N2P")
  }

  return(result)
}