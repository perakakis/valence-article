# AR.R: Compute Autoregressive coefficient per subject using mixed-effects model
# Translation of MATLAB statistics/AR.m
#
# Function: AR
# Inputs:
#   data: Data.table with columns PID, day, beep, and the variables specified in vars
#   vars: Character vector of column names (e.g., c("PosA", "NegA"))
# Outputs:
#   Data.table with columns PID, AR_<var1>, AR_<var2>, ...
#
# Algorithm:
#   1. Person-mean center the data
#   2. Fit mixed model: X(t+1) ~ X(t) + (1 + X(t) | PID)
#   3. Extract per-subject AR from fixed + random effects

AR <- function(data, vars) {
    library(data.table)
    library(lme4)

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

    # For each variable, compute AR per subject using mixed model
    for (v in vars) {
        col_name <- paste0("AR_", v)
        result[, (col_name) := NA_real_]

        # Person-mean center the variable (matching MATLAB line 18-19)
        dt[, centered := get(v) - mean(get(v), na.rm = TRUE), by = PID]

        # Create lagged pairs: X(t) and X(t+1) for consecutive observations
        # Conditions: same PID, same day, beep+1
        dt[, row_id := .I]
        dt[, next_centered := shift(centered, type = "lead"), by = .(PID, day)]
        dt[, next_beep := shift(beep, type = "lead"), by = .(PID, day)]

        # Filter to valid consecutive pairs
        pairs <- dt[!is.na(centered) & !is.na(next_centered) &
                    next_beep == beep + 1,
                    .(PID, X_t = centered, X_t1 = next_centered)]

        if (nrow(pairs) < 10) {
            warning("Not enough consecutive pairs for AR estimation of ", v)
            next
        }

        # Fit mixed model: X(t+1) ~ X(t) + (1 + X(t) | PID)
        # This matches MATLAB: mdl = fitlmematrix([1, X(t)], X(t+1), [1, X(t)], subj)
        tryCatch({
            mdl <- lmer(X_t1 ~ X_t + (1 + X_t | PID), data = pairs,
                       control = lmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 10000)))

            # Extract per-subject AR coefficients
            # Fixed effect + random effect for slope
            fixed_slope <- fixef(mdl)["X_t"]
            random_slopes <- ranef(mdl)$PID$X_t

            # Get AR for each PID that was in the model
            model_pids <- rownames(ranef(mdl)$PID)
            for (i in seq_along(model_pids)) {
                pid <- as.numeric(model_pids[i])
                ar_value <- fixed_slope + random_slopes[i]
                result[PID == pid, (col_name) := ar_value]
            }
        }, error = function(e) {
            warning("AR model failed for ", v, ": ", e$message)
        })

        # Clean up temporary columns
        dt[, c("centered", "row_id", "next_centered", "next_beep") := NULL]
    }

    return(result)
}