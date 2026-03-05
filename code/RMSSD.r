# RMSSD.R: Compute Root Mean Square of Successive Differences per subject
# Translation of MATLAB statistics/MSSD.m and statistics/RMSSD.m
#
# Function: RMSSD
# Inputs:
#   data: Data.table with columns PID, day, beep, and the variables specified in vars
#   vars: Character vector of column names (e.g., c("PosA", "NegA"))
# Outputs:
#   Data.table with columns PID, RMSSD_<var1>, RMSSD_<var2>, ...
#
# Algorithm:
#   For consecutive observations (t, t+1) within same subject, day, and beep+1:
#   MSSD = mean((X(t+1) - X(t))^2)
#   RMSSD = sqrt(MSSD)

RMSSD <- function(data, vars) {
    library(data.table)

    # Validate inputs
    if (!"PID" %in% names(data)) stop("PID column missing in data")
    if (!"day" %in% names(data)) stop("day column missing in data")
    if (!"beep" %in% names(data)) stop("beep column missing in data")
    for (v in vars) {
        if (!v %in% names(data)) stop("Missing variable in data: ", v)
    }

    # Copy data (assumes data is already sorted by PID, day, beep)
    dt <- copy(data)

    # Initialize result
    pids <- unique(dt$PID)
    result <- data.table(PID = pids)

    # For each variable, compute RMSSD per subject
    for (v in vars) {
        col_name <- paste0("RMSSD_", v)
        result[, (col_name) := NA_real_]

        for (s in pids) {
            # Get subject data
            subj_dt <- dt[PID == s]
            n <- nrow(subj_dt)

            if (n < 2) {
                next  # Not enough data
            }

            # Find consecutive pairs: same day, beep+1
            # Compare row i with row i+1
            sq_diffs <- numeric(0)

            for (i in 1:(n - 1)) {
                # Check conditions:
                # 1. Same day
                # 2. Consecutive beep (beep[i+1] == beep[i] + 1)
                # 3. Both values non-NA
                if (subj_dt$day[i] == subj_dt$day[i + 1] &&
                    subj_dt$beep[i + 1] == subj_dt$beep[i] + 1 &&
                    !is.na(subj_dt[[v]][i]) &&
                    !is.na(subj_dt[[v]][i + 1])) {

                    diff_sq <- (subj_dt[[v]][i + 1] - subj_dt[[v]][i])^2
                    sq_diffs <- c(sq_diffs, diff_sq)
                }
            }

            # Compute RMSSD if we have valid pairs
            if (length(sq_diffs) > 0) {
                mssd <- mean(sq_diffs)
                rmssd <- sqrt(mssd)
                result[PID == s, (col_name) := rmssd]
            }
        }
    }

    return(result)
}
