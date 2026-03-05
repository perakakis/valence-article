# Load necessary libraries
library(dplyr)
library(psych)
options(warn = -1)  # Suppress warnings for cleaner output

# Load external auxiliary functions
source("./code/auxFunctions.R")

# Load and process data for Study 1 (Spanish EMA data)
# df1: Spanish study ####
df1_ema <- read.csv(file = "./data/spanish_EMA.csv")  # Load EMA data
df1_ema <- rename(df1_ema, day = study_day, beep = occasion)  # Rename columns for clarity

# Replace rows where valence, arousal, and energy remained at initial values (likely invalid data)
ind <- which(df1_ema$valence == 0 & df1_ema$energetic_arousal == 50 & df1_ema$tense_arousal == 50)
df1_ema$valence[ind] = NA  # Mark those rows as NA for valence

# Estimate adherence for each participant
PIDs <- unique(df1_ema$participant)  # Get unique participant IDs

for (i in seq(1, length(PIDs))) {
  y <- subset(df1_ema, participant == PIDs[i])  # Subset data by participant
  y <- y$valence  # Select valence variable
  # Calculate adherence as the proportion of non-missing valence entries
  df1_ema$adh[df1_ema$participant == PIDs[i]] <- 1 - length(which(is.na(y))) / length(y)
}

# Load and process survey data for Study 1
df1_survey <- read.csv(file = "./data/spanish_survey.csv")  # Load survey data

# Compute total scores for psychological measures
df1_survey$GAD <- rowSums(df1_survey[, c("GAD_1", "GAD_2", "GAD_3", "GAD_4", "GAD_5", "GAD_6", "GAD_7")])
df1_survey$PHQ <- rowSums(df1_survey[, c("PHQ_1", "PHQ_2", "PHQ_3", "PHQ_4", "PHQ_5", "PHQ_6", "PHQ_7", "PHQ_8", "PHQ_9")])
df1_survey$BRS <- rowMeans(cbind(df1_survey$BRS_1, df1_survey$BRS_3, df1_survey$BRS_5, 6 - df1_survey$BRS_2, 6 - df1_survey$BRS_4, 6 - df1_survey$BRS_6))
df1_survey$FS <- rowSums(df1_survey[, paste0("FS_", 1:8)])
df1_survey$AAQ <- rowSums(df1_survey[, paste0("AAQ_", 1:7)])

# Rename columns for consistency and clarity
df1_survey <- rename(df1_survey, SWLS = GLS, sex = gender)

# Merge EMA and survey data for Study 1
df1 <- merge(df1_ema, df1_survey, by = "participant")
df1$Study <- 1  # Add a study identifier

# Create a participant number vector
df1 <- transform(df1, PID = match(participant, unique(participant)))

# Select relevant columns and sort by participant, day, and beep
df1 <- df1[c("Study", "participant", "PID", "age", "sex", "adh", "day", "beep", "valence", "PHQ", "GAD", "AAQ", "FS", "BRS", "SWLS")]
df1 <- arrange(df1, PID, day, beep)

# Load and process data for Study 2 (German EMA data)
# df2: German study ####
df2_ema <- read.csv(file = "./data/german_EMA.csv")
df2_ema <- rename(df2_ema, PID = VID, valence = Core.Affect.Valence, arousal = Core.Affect.Arousal, day = Day, beep = Session)

# Replace rows where valence and arousal remained at initial values
ind <- which(df2_ema$valence == 50 & df2_ema$arousal == 50)
df2_ema$valence[ind] = NA  # Mark those rows as NA

# Replace invalid valence values (-1) with NA
ind <- which(df2_ema$valence == -1)
df2_ema$valence[ind] = NA

# Estimate adherence for each participant in Study 2
PIDs <- unique(df2_ema$PID)
for (i in seq(1, length(PIDs))) {
  y <- subset(df2_ema, PID == PIDs[i])
  y <- y$valence  # Select valence variable
  df2_ema$adh[df2_ema$PID == PIDs[i]] <- 1 - length(which(is.na(y))) / length(y)
}

# Adjust valence values by subtracting 50 for scaling purposes
df2_ema$valence <- df2_ema$valence - 50

# Load and merge survey data for Study 2
df2_survey <- read.csv(file = "./data/00 german_survey.csv")
df2_survey <- rename(df2_survey, PID = VID)
df2 <- merge(df2_ema, df2_survey, by = "PID")
df2$Study <- 2  # Add a study identifier

# Create unique participant IDs, continuing from Study 1
df2 <- transform(df2, id = match(PID, unique(PID)))
df2$id <- df2$id + tail(df1$PID, n = 1)

# Preserve the original participant variable
df2$participant <- df2$PID

# Select relevant columns and rename for consistency
df2 <- df2[c("Study", "participant", "id", "adh", "day", "beep", "Gender.x", "Age.x", "valence", "DASS.Depression", "DASS.Anxiety", "AAQ", "BriefResilience", "SatisfactLife")]
df2 <- rename(df2, PID = id, DASSd = DASS.Depression, DASSa = DASS.Anxiety, SWLS = SatisfactLife, sex = Gender.x, age = Age.x, BRS = BriefResilience)

# Combine data from both studies
# df1
N1 <- length(unique(df1$PID))  # Number of participants in Study 1
df1 <- subset(df1, df1$adh > 0.75)  # Apply adherence criterion (>75%)
# Apply additional exclusion criteria using auxiliary function
exclusion_list <- identify_exclusions(df1, "PID", "valence", 20)
df1 <- subset(df1, !PID %in% exclusion_list)
N2 <- length(unique(df1$PID))  # Number of participants after exclusions

# df2
N1 <- length(unique(df2$PID))  # Number of participants in Study 2
df2 <- subset(df2, df2$adh > 0.75)  # Apply adherence criterion (>75%)
# Apply additional exclusion criteria
exclusion_list <- identify_exclusions(df2, "PID", "valence", 20)
df2 <- subset(df2, !PID %in% exclusion_list)
N2 <- length(unique(df2$PID))  # Number of participants after exclusions

# Combine datasets from both studies
dfs <- plyr::rbind.fill(df1, df2)  # Combine data
dfs <- dfs[c("Study", "participant", "PID", "age", "sex", "adh", "day", "beep", "valence", "GAD", "PHQ", "AAQ", "FS", "BRS", "SWLS", "DASSd", "DASSa")]

# Final count of participants across both studies
totalN <- length(unique(dfs$PID))  # Total number of participants after exclusions