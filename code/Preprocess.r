# Preprocess.R: Load, preprocess, compute stats for EMA datasets
# Outputs SAll, RAll, CAll,COVDataset, infoDatasets
library(data.table)

# Create output directories if needed
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("plots")) dir.create("plots")

# Load filter participants script
source("code/FilterParticipants.R")

# Load statistics scripts
source("code/MEAN_SD.R")
source("code/ASR.R")
source("code/RMSSD.r")
source("code/AR.r")

DataStrAll <- c("00 german", "01 postcovid1", "02 postcovid2")

# Select outcomes
outcome_groups <- list(
    Depression = c("CESD", "QIDS", "PHQ", "BDI", "DASSd"),
    Anxiety = c("GAD", "DASSa"),
    Other = c("SWL", "BRS", "AAQ", "FS", "SWL1")
)

# Parameters
minCompliance <- 0.75
center_value <- 50 # use to center valence if necessary

# Column mappings
col_names <- list(
    PID = "PID",
    day = "UNIT",
    beep = "BEEP",
    occasion = "OCCASION",
    PosA = "PA",
    NegA = "NA",
    valence = "Happy"
)

# Initialize lists
SAll <- vector("list", length(DataStrAll))
RAll <- vector("list", length(DataStrAll))
ValenceDataAll <- vector("list", length(DataStrAll))

# Loop over datasets
for (i in seq_along(DataStrAll)) {
    dataStr <- DataStrAll[i]
    cat(sprintf("Processing %s (%d/%d)\n", dataStr, i, length(DataStrAll)))

    # Load main data
    file <- paste0("data/", dataStr, ".csv")
    if (!file.exists(file)) stop("File not found: ", file)
    data <- fread(file,
        sep = "auto", header = TRUE, na.strings = c(""),
        colClasses = list(character = "PID"), check.names = FALSE
    )

    # Check all columns exist
    required_cols <- unlist(col_names)
    if (!all(required_cols %in% names(data))) {
        stop("Missing columns in ", file, ": ", paste(setdiff(required_cols, names(data)), collapse = ", "))
    }

    # Ensure numeric columns
    for (col in required_cols) {
        if (is.character(data[[col]]) && col != "PID") {
            data[, (col) := as.numeric(get(col))]
        }
    }

    # Extract relevant columns and rename
    data <- data[, ..required_cols]
    setnames(data, names(col_names))

    # Load regression outcomes
    file <- paste0("data/", dataStr, " reg.csv")
    if (file.exists(file)) {
        Rfull <- fread(file,
            sep = "auto", header = TRUE, na.strings = c("", "NA"),
            colClasses = list(character = "PID")
        )
    } else {
        Rfull <- NULL
    }

    # Preprocess dataset
    preprocess_result <- FilterParticipants(
        data = data,
        Rfull = Rfull,
        # Cfull = Cfull,
        minCompliance = minCompliance
        # filterConsecMissing = filterConsecMissing
    )

    data <- preprocess_result$dataNew
    Rfull <- preprocess_result$RfullNew

    # Compute stats
    S_stat <- MEAN_SD(data, c("PosA", "NegA"))
    S_asr <- ASR(data, "valence", center_value = center_value)
    S_rmssd <- RMSSD(data, c("PosA", "NegA"))
    S_ar <- AR(data, c("PosA", "NegA"))

    # Merge all stats
    S_merged <- merge(S_stat, S_asr, by = "PID")
    S_merged <- merge(S_merged, S_rmssd, by = "PID")
    S_merged <- merge(S_merged, S_ar, by = "PID")
    SAll[[i]] <- S_merged
    RAll[i] <- list(Rfull)

    # Extract valence timeseries data for export
    ValenceDataAll[[i]] <- data[, .(original_PID = PID, day = day, beep = beep, valence = valence, dataset = dataStr)]

    cat("Processing complete.\n")
}
# Build long-format table with outcome groups
outcome_type_lookup <- data.table(
    outcome_name = unlist(outcome_groups, use.names = FALSE),
    outcome_type = rep(names(outcome_groups), times = lengths(outcome_groups))
)

long_data_list <- vector("list", length(DataStrAll))
for (i in seq_along(DataStrAll)) {
    dataset_name <- DataStrAll[i]
    S <- SAll[[i]]
    R <- RAll[[i]]

    # Regression outcomes
    if (!is.null(S) && !is.null(R)) {
        merged <- merge(R, S, by = "PID", all = FALSE)
        outcome_names <- setdiff(names(R), "PID")
        melted <- suppressWarnings(melt(
            merged,
            id.vars = setdiff(names(merged), outcome_names),
            measure.vars = outcome_names,
            variable.name = "outcome_name",
            value.name = "outcome_value"
        ))
        melted <- as.data.table(melted)
        melted[, dataset := dataset_name]
        # Add outcome_type using the lookup table
        melted <- merge(melted, outcome_type_lookup, by = "outcome_name", all.x = TRUE)
        unknown <- unique(melted[is.na(outcome_type), outcome_name])
        if (length(unknown) > 0) {
            warning(sprintf(
                "Dataset '%s': outcome(s) not in outcome_groups and will be ignored: %s",
                dataset_name, paste(unknown, collapse = ", ")
            ))
        }
        long_data_list[[i]] <- melted
    }
}
long_data <- rbindlist(long_data_list, fill = TRUE)

# Clean up long_data - convert outcome_value to numeric if it's character
if (is.character(long_data$outcome_value)) {
    long_data[, outcome_value := as.numeric(outcome_value)]
}

# Reorder columns
if (nrow(long_data) > 0) {
    setcolorder(long_data, c(
        "dataset", "PID", "outcome_name",
        "outcome_type", "outcome_value"
    ))
}

# Combine RAll and SAll into a single long format CSV
# Bind all regression outcomes with dataset names
RAll_combined <- rbindlist(RAll, idcol = "dataset_id", fill = TRUE)
if (nrow(RAll_combined) > 0) {
    RAll_combined[, dataset := DataStrAll[dataset_id]]
    RAll_combined[, dataset_id := NULL]
}

# Bind all statistics with dataset names
SAll_combined <- rbindlist(SAll, idcol = "dataset_id", fill = TRUE)
if (nrow(SAll_combined) > 0) {
    SAll_combined[, dataset := DataStrAll[dataset_id]]
    SAll_combined[, dataset_id := NULL]
}

# Merge RAll and SAll by dataset and PID
if (nrow(RAll_combined) > 0 && nrow(SAll_combined) > 0) {
    combined_data <- merge(RAll_combined, SAll_combined, by = c("dataset", "PID"), all = TRUE)
} else if (nrow(SAll_combined) > 0) {
    combined_data <- SAll_combined
} else if (nrow(RAll_combined) > 0) {
    combined_data <- RAll_combined
} else {
    combined_data <- data.table()
}

# Save combined data
if (nrow(combined_data) > 0) {
    write.csv(combined_data, "results/AllData_wide.csv", row.names = FALSE)
    cat(sprintf("Combined RAll and SAll saved to: results/AllData_wide.csv (%d rows)\n", nrow(combined_data)))
} else {
    cat("Warning: No data to save in combined RAll and SAll file.\n")
}

# Create valence timeseries CSV with PID index and time counter
cat("\nCreating valence timeseries export...\n")

# Combine all valence data
valence_combined <- rbindlist(ValenceDataAll, fill = TRUE)

if (nrow(valence_combined) > 0) {
    valence_combined[, composite_key := paste(dataset, original_PID, sep = "_")]
    unique_keys <- unique(valence_combined$composite_key)
    pid_mapping <- data.table(composite_key = unique_keys, PID = seq_along(unique_keys))
    valence_combined <- merge(valence_combined, pid_mapping, by = "composite_key")
    setorder(valence_combined, PID, day, beep)
    valence_combined[, time := seq_len(.N), by = PID]
    valence_export <- valence_combined[, .(PID, time, day, valence, dataset)]

    # Export to CSV
    write.csv(valence_export, "results/valence_timeseries.csv", row.names = FALSE)
    cat(sprintf(
        "Valence timeseries saved to: results/valence_timeseries.csv (%d rows, %d unique PIDs)\n",
        nrow(valence_export), length(unique(valence_export$PID))
    ))
} else {
    cat("Warning: No valence data to export.\n")
}