# Load libraries, define parameters ----
library(dplyr)
library(earlywarnings)
library(photobiology)

# Clear workspace
rm(list = ls())

# Define parameters
pl <- FALSE # Set to TRUE to generate plots
boundsAll <- c(0, 100)
emotion_variable <- c("Happy", "happy") # The variable to analyze
minCompliance <- 0.5 # Minimum compliance for inclusion in the analysis

# List to store all stats and outcomes
allData <- list()

# All datasets
# dataStrAll <- c(
#     "MDD GOTLIB", "KATHLEEN STRESSKLINIEK", "MTURK DAILY DIARY",
#     "Cogito", "Laura ESM 2016", "Laura ESM 2014", "MARLIES BPD",
#     "PETEMADDY", "CLINICAL ESM", "ELISE ESM14", "LONGITUDINAL W1",
#     "LONGITUDINAL W2", "LONGITUDINAL W3", "MDD BPD TRULL", "JULIAN EGON"
# )

# Exclude MDD GOTLIB and MDD BPD TRULL because they do not have a reg file.
dataStrAll <- c(
    "KATHLEEN STRESSKLINIEK", "MTURK DAILY DIARY",
    "Cogito", "Laura ESM 2016", "Laura ESM 2014", "MARLIES BPD",
    "PETEMADDY", "CLINICAL ESM", "ELISE ESM14", "LONGITUDINAL W1",
    "LONGITUDINAL W2", "LONGITUDINAL W3", "JULIAN EGON"
)

dataStrAll <- c(
    "KATHLEEN STRESSKLINIEK", 
    "Cogito", "Laura ESM 2016", "Laura ESM 2014", "MARLIES BPD",
    "PETEMADDY", "CLINICAL ESM", "ELISE ESM14", "LONGITUDINAL W1",
    "LONGITUDINAL W2", "LONGITUDINAL W3", "JULIAN EGON"
)

# dataStrAll <- c(
#     "PETEMADDY", "ELISE ESM14"
# )

# Mapping of dataset codenames to paper references
codenames2papers <- data.frame(
  codename = c(
    "MDD GOTLIB", "KATHLEEN STRESSKLINIEK", "MTURK DAILY DIARY",
    "Cogito", "Laura ESM 2016", "Laura ESM 2014", "MARLIES BPD",
    "PETEMADDY", "CLINICAL ESM", "ELISE ESM14", "LONGITUDINAL W1",
    "LONGITUDINAL W2", "LONGITUDINAL W3", "MDD BPD TRULL", "JULIAN EGON"
  ),
  paper = c(
    "Thompson 2012", "Van der Gucht 2017", "Dejonckheere 2017",
    "Schmiedek 2010", "Sels 2018", "Sels 2017", "Houben 2016",
    "Koval 2013", "Heininga 2019", "Dejonckheere 2018", "Pe 2016 Wave 1",
    "Pe 2016 Wave 2", "Pe 2016 Wave 3", "Trull 2008", "Provenzano in prep"
  ),
  stringsAsFactors = FALSE
)

# Dataset loop ----
for (dataStr in dataStrAll) {
    ## Load the dataset ----
    df_esm <- read.csv(paste0("./data/", dataStr, ".csv"), sep = ";")
    df_reg <- read.csv(paste0("./data/", dataStr, " reg.csv"), sep = ";")
    bounds <- read.csv(paste0("./data/", dataStr, " bound.csv"), sep = ";")

    # Check if the emotion variable exists; if not, skip this dataset
    if (!any(emotion_variable %in% colnames(df_esm))) {
        cat("Skipping dataset", dataStr, "- emotion variable not found\n")
        next
    }

    # Rescale emotion variables (from 4:end) to 0-100 ----
    df_esm[, 4:ncol(df_esm)] <- (df_esm[, 4:ncol(df_esm)] - bounds$LOWER) / (bounds$UPPER - bounds$LOWER) * (boundsAll[2] - boundsAll[1]) + boundsAll[1]

    ## Subject loop ----
    dataset_stats <- data.frame()
    PIDs <- unique(df_esm$PID)

    for (i in 1:length(PIDs)) {
        # Select the data for the current subject ----
        data <- df_esm[df_esm$PID == PIDs[i], ]

        # Sort data by OCCASION
        data <- data[order(data$OCCASION), ]

        # Get emotion variable
        y <- data[[emotion_variable[emotion_variable %in% colnames(data)]]]
        y <- y - 50 # Transform to bipolar

        # Check minCompliance
        valid_responses <- sum(!is.na(y)) / length(y)
        if (valid_responses < minCompliance) {
            cat("Skipping subject", PIDs[i], "in dataset", dataStr, "- insufficient valid responses\n")
            next
        }

        # Compute statistics ----
        ## Mean and SD ----
        mPA <- mean(data$PA, na.rm = TRUE)
        mNA <- mean(data$NA., na.rm = TRUE)
        sdPA <- sd(data$PA, na.rm = TRUE)
        sdNA <- sd(data$NA., na.rm = TRUE)

        ## Bistability classification ----

        # Create a valence variable without NAs ----
        y_clean <- y[!is.na(y)]
        basins <- earlywarnings::livpotential_ews(y_clean)

        # Identify local minima (valleys) in the potential function
        lm <- photobiology::get_valleys(
            basins$grid.points, basins$pot,
            ignore_threshold = -0.2,
            strict = TRUE,
            span = 3
        )

        # Extract indices of local minima (valleys)
        basinsind <- which(basins$grid.points %in% lm$x)
        basinlocs <- lm$x # Coordinates of basin locations

        # Classify the stability landscape
        if (length(basinlocs) == 1) {
            if (basinlocs >= 0) {
                type <- "positive monostable" # One local minimum, positive valence
            } else {
                type <- "negative monostable" # One local minimum, negative valence
            }
        } else if (length(basinlocs) > 1) {
            if (all(basinlocs >= 0)) {
                type <- "positive multistable" # Multiple local minima, all positive
            } else if (all(basinlocs <= 0)) {
                type <- "negative multistable" # Multiple local minima, all negative
            } else if (any(basinlocs > 0) & any(basinlocs < 0)) {
                type <- "bistable" # Mixed positive and negative local minima
            }
        } else {
            type <- "undefined" # No local minima found
        }

        # Classify bistability (1 = bistable, 0 = monostable)
        if (any(basinlocs > 0) && any(basinlocs < 0)) {
            bistable <- 1 # Bistable
        } else {
            bistable <- 0 # Monostable
        }

        ## Affect shift metrics ----
        ### ASR metrics ----
        p2n <- which(diff(sign(y)) == -2) # Positive-to-negative shifts
        n2p <- which(diff(sign(y)) == 2) # Negative-to-positive shifts

        positive_obs <- length(which(y > 0))
        if (positive_obs > 0) {
            if (length(p2n) > 0) {
                P2N <- length(p2n) / positive_obs
                mP2N <- mean(y[p2n] + abs(y[p2n + 1]))
                sdP2N <- sd(y[p2n] + abs(y[p2n + 1]))
            } else {
                P2N <- 0
                mP2N <- NA
                sdP2N <- NA
            }
        } else {
            P2N <- 1 #NA# 1 # Change to NA for accurate interpretation 
            mP2N <- NA
            sdP2N <- NA
        }

        negative_obs <- length(which(y < 0))
        if (negative_obs > 0) {
            if (length(n2p) > 0) {
                N2P <- length(n2p) / negative_obs
                mN2P <- mean(abs(y[n2p]) + y[n2p + 1])
                sdN2P <- sd(abs(y[n2p]) + y[n2p + 1])
            } else {
                N2P <- 0
                mN2P <- NA
                sdN2P <- NA
            }
        } else {
            N2P <- 0 #NA #0 # Change to NA for accurate interpretation
            mN2P <- NA
            sdN2P <- NA
        }

        # Residence Time metrics ----
        PAt <- rle(sign(y))$lengths[(rle(sign(y))$values == 1)] # Positive affect residence time
        NAt <- rle(sign(y))$lengths[(rle(sign(y))$values == -1)] # Negative affect residence time

        # Positive residence time statistics
        if (length(PAt) > 0 & !all(is.na(PAt))) {
            mPRT <- mean(PAt, na.rm = TRUE)
            sdPRT <- sd(PAt, na.rm = TRUE)
        }

        # Negative residence time statistics
        if (length(NAt) > 0 & !all(is.na(NAt))) {
            mNRT <- mean(NAt, na.rm = TRUE)
            sdNRT <- sd(NAt, na.rm = TRUE)
        }

        # Subject results
        subject_stats <- data.frame(
            Dataset = dataStr,
            PID = PIDs[i],
            mPA = mPA,
            mNA = mNA,
            sdPA = sdPA,
            sdNA = sdNA,
            P2N = P2N,
            N2P = N2P,
            mP2N = mP2N,
            sdP2N = sdP2N,
            mN2P = mN2P,
            sdN2P = sdN2P,
            mPRT = mPRT,
            sdPRT = sdPRT,
            mNRT = mNRT,
            sdNRT = sdNRT,
            Bistable = bistable,
            Type = type
        )

        # Merge with reg data
        reg_vars <- df_reg[df_reg$PID == PIDs[i], ]
        reg_vars <- reg_vars[, names(reg_vars) != "PID", drop = FALSE]
        subject_stats <- cbind(subject_stats, reg_vars)

        # Update the dataset statistics
        dataset_stats <- rbind(dataset_stats, subject_stats)

        # Generate plots if the plot flag is set to TRUE
        if (pl) {
            # Create folder for dataset-specific figures if it doesn't exist
            dataset_fig_dir <- file.path("./figures", dataStr, "subjects")
            if (!dir.exists(dataset_fig_dir)) {
                dir.create(dataset_fig_dir, recursive = TRUE)                
            }

            pdf(paste0(dataset_fig_dir, "/", PIDs[i], ".pdf"), width = 8, height = 8)
            par(mfrow = c(2, 1), mar = c(5, 3, 3, 3))

            # Plot valence over time
            plot(y, type = "l", ylim = c(-50, 50), main = paste0("PID: ", PIDs[i]), xlab = "", ylab = "")

            # Plot histograms and basins
            h <- hist(y, breaks = 50, plot = FALSE)
            # plot(h, main = paste0("type = ", type), xlim = c(-50, 50), xlab = "")
            plot(h, xlim = c(-50, 50), xlab = "")

            # par(new = TRUE)
            # plot(basins$grid.points, basins$pot, xlim = c(-50, 50), axes = FALSE, xlab = "", ylab = "", type = "l")
            # axis(side = 4, at = pretty(range(basins$pot)))
            # par(new = TRUE)
            # plot(basins$grid.points[basinsind], basins$pot[basinsind], col = 2, xlim = c(-50, 50), axes = FALSE, xlab = "", ylab = "", ylim = range(basins$pot))

            # Close the PDF file for this participant
            dev.off()
        }
    }
    # Update all stats
    allData[[dataStr]] <- dataset_stats
}

combined_results <- dplyr::bind_rows(allData)
write.csv(combined_results, "./results/allData.csv", row.names = FALSE)

# Save the allData to a file
saveRDS(allData, "./results/allData.rds")