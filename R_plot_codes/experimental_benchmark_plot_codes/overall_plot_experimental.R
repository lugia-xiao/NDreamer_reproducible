# Load necessary library for string manipulation
library(dplyr)

# Define the function
normalize_and_aggregate <- function(dataframe) {
  # Rename columns to correct issues
  colnames(dataframe) <- colnames(dataframe) %>% 
    gsub("X1.cLISI", "cLISI", .) %>% 
    gsub("kBET.Accept.Rate", "kBET", .) %>% 
    gsub("True Pos Rate", "True.Pos.Rate", .)
  
  # Set environment language
  Sys.setenv(LANG = "en")
  
  # Normalize the float columns using 01 normalization
  numeric_cols <- setdiff(names(dataframe), "Method")
  dataframe[numeric_cols] <- lapply(dataframe[numeric_cols], function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  
  # Compute mean scores for each group
  group1 <- c("bLISI", "ASW_batch", "kBET", "True.Pos.Rate")
  group2 <- c("ASW_label", "ARI", "NMI", "cLISI")
  
  # Ensure columns exist in the dataframe
  if (!all(group1 %in% names(dataframe)) || !all(group2 %in% names(dataframe))) {
    stop("One or more columns in the specified groups are missing in the dataframe.")
  }
  
  dataframe$Score1 <- rowMeans(dataframe[group1], na.rm = TRUE)
  dataframe$Score2 <- rowMeans(dataframe[group2], na.rm = TRUE)
  
  # Return the dataframe with scores
  return(dataframe[, c("Method", "Score1", "Score2")])
}

# Set working directory
setwd("D:/project/ndreamer_plot/")

# Read the data and apply the function
paths <- c("./experimental/PBMC.csv", "./experimental/virus.csv", "./experimental/PBMC_yale.csv")

for (pathi in paths) {
  # Read the data
  df <- read.csv(pathi)
  
  # Apply the normalization and aggregation function
  result <- normalize_and_aggregate(df)
  
  # Print the result
  print(paste("Results for file:", pathi))
  print(result)
}

# Create the three dataframes
data_pbmc <- data.frame(
  Method = c("NDreamer", "CINEMA-OT", "Mixscape", "scGen", "scCAPE"),
  Score1 = c(1.00000000, 0.31980917, 0.10562393, 0.06953034, 0.33832554),
  Score2 = c(0.8081643, 0.1371153, 0.4672925, 0.3579749, 0.2130101)
)

data_virus <- data.frame(
  Method = c("NDreamer", "CINEMA-OT", "Mixscape", "scGen", "scCAPE"),
  Score1 = c(1.00000000, 0.47799963, 0.02892542, 0.58730898, 0.07157991),
  Score2 = c(0.6999977, 0.2214840, 0.9275032, 0.4649229, 0.1172412)
)

data_pbmc_yale <- data.frame(
  Method = c("NDreamer", "CINEMA-OT", "Mixscape", "scGen", "scCAPE"),
  Score1 = c(0.99810671, 0.65123946, 0.07818016, 0.28712101, 0.65426330),
  Score2 = c(0.7431584, 0.4123340, 0.8661154, 0.3158624, 0.4420572)
)

# Combine the dataframes and calculate the mean scores
all_data <- rbind(data_pbmc, data_virus, data_pbmc_yale)
averaged_scores <- aggregate(. ~ Method, data = all_data, FUN = mean)

# Scatter plot with customizations
library(ggplot2)

ggplot(averaged_scores, aes(x = Score1, y = Score2, label = Method, color = Method)) +
  geom_point(size = 3) +  # Increase point size
  geom_text(vjust = -0.5, hjust = 0.5, fontface = "bold") +  # Add labels
  ggtitle("Condition Mixing vs. Biological Conservation") +
  xlab("Average Score (Condition Mixing)") +
  ylab("Average Score (Biological Conservation)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"  # Remove legend if not needed
  )