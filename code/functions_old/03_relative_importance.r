# Load external scripts for data extraction, auxiliary functions, and custom color schemes
source("./code/01_preprocess_extract.R")  
source("./code/auxFunctions.R")
source("./code/predictor_colors.R")

# Suppress warnings for cleaner output
options(warn = -1)

# Load necessary libraries
library(relaimpo)   # For calculating relative importance
library(reshape2)   # For reshaping data
library(ggsci)      # For color palettes
library(ggplot2)    # For plotting

# Set seed for reproducibility
set.seed(1)

# Subset the data for Study 1 and aggregate to get only the first row per participant
df <- aggregate(df, list(df$PID), data = df, FUN = head, 1)

# Select the relevant metrics for relative importance analysis
metrics <- df[c("mPA", "sdPA", "mNA", "sdNA", "P2N_ASR","N2P_ASR")]  # Positive and negative affect, shifts, and standard deviations

# Define outcomes of interest
outcomes <- c("CESD", "SWL")

# Initialize an empty dataframe to store the relative importance results for Study 1
lres <- data.frame(matrix(ncol = length(metrics), nrow = length(outcomes)))  # Create empty dataframe
colnames(lres) <- colnames(metrics)  # Set column names to metric names
rownames(lres) <- outcomes  # Set row names to outcome names

# Perform relative importance analysis for each outcome in Study 1
for (i in 1:length(outcomes)) {
  outcome <- outcomes[i]  # Define the outcome variable (e.g., GAD, PHQ)
  res <- relativeImportance(df[[outcome]], metrics)  # Perform relative importance analysis using custom function
  lres[i, names(res$lmg)] <- res$lmg  # Store the relative importance (lmg)
}

data <- lres
colnames(data)[5:6] <- c("P2N-ASR","N2P-ASR")  # Rename columns for clarity

# Add the 'outcome' as a separate column for the plot
data$outcome <- rownames(data)

# Reshape data for plotting
data_long <- melt(data, id.vars = "outcome", variable.name = "predictor", value.name = "importance")

# Define x-axis labels for the plot
xlabels <- c("CESD", "SWL")

# Create relative importance plot
p <- ggplot(data_long, aes(x = outcome, y = importance, fill = predictor)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +  # Plot importance as bars
  scale_fill_manual(values = predictor_colors, name = "Predictor") +  # Use custom color palette for predictors
  scale_x_discrete(limits = unique(data_long$outcome), labels = xlabels) +  # Customize x-axis labels
  theme_minimal() +  # Use a minimal theme for a clean plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(title = "Relative importance",
       x = "Well-being Indicator",
       y = "Relative Importance") +  # Add titles and labels to the plot
  theme(axis.title = element_text(size = 24),  # Customize axis title size
        axis.text = element_text(size = 20, colour = "black"),  # Customize axis text size and color
        plot.title = element_text(size = 24),  # Customize plot title size
        legend.text = element_text(size = 16),  # Customize legend text size
        legend.title = element_text(size = 16))  # Customize legend title size

# Optional: Save the plot as a PDF
# ggsave("./figures/relImportance.pdf", p, width = 10, height = 6, units = "in")