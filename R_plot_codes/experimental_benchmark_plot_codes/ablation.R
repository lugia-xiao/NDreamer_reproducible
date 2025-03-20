# Set environment language
Sys.setenv(LANG = "en")

# Set working directory
setwd("D:/project/ndreamer_plot/")
library(ggplot2)

# Read the data
path <- c("./experimental/ablation.csv")
df <- read.csv(path)

# Rename columns
colnames(df) <- gsub("^X1\\.cLISI$", "cLISI", colnames(df))
colnames(df) <- gsub("^kBET\\.Accept\\.Rate$", "kBET", colnames(df))

# Define metric groups
batch_correction_metrics <- c("bLISI", "ASW_batch", "kBET", "True.Pos.Rate")
biological_conservation_metrics <- c("ASW_label", "ARI", "NMI", "cLISI")

# Plot batch correction metrics
batch_correction_data <- df[, c("Model", batch_correction_metrics)]
batch_correction_long <- reshape2::melt(batch_correction_data, id.vars = "Model")

ggplot(batch_correction_long, aes(x = Model, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Condition-mixing Metrics",
       x = "Model Configuration",
       y = "Metric Value",
       fill = "Metrics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot biological conservation metrics
biological_conservation_data <- df[, c("Model", biological_conservation_metrics)]
biological_conservation_long <- reshape2::melt(biological_conservation_data, id.vars = "Model")

ggplot(biological_conservation_long, aes(x = Model, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Biological Conservation Metrics",
       x = "Model Configuration",
       fill = "Metrics") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
