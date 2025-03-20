# Set environment language
Sys.setenv(LANG = "en")

# Set working directory
setwd("D:/project/ndreamer_plot/")

# Load necessary libraries
library(dplyr)
library(fmsb)
library(gridExtra)
library(gt)
library(ggplot2)

generate_radar_plot <- function(df,title) {
  # Rename columns
  colnames(df) <- gsub("^X1\\.cLISI$", "cLISI", colnames(df))
  colnames(df) <- gsub("kBET\\.Accept\\.Rate", "kBET", colnames(df))
  
  # Group definitions
  groups <- list(
    "EM-batch condition mix" = c("bLISI_EM", "ASW_batch_EM", "kBET_EM", "True.Pos.Rate_EM"),
    "Denoised-batch mixing" = c("bLISI_denoised", "ASW_batch_denoised", "kBET_denoised", "True.Pos.Rate_denoised"),
    "Signal preservation" = c("trace_similarity"),
    "Data distortion" = c("F1", "Accuracy", "mean_r")
  )
  
  # Step 1: Calculate the mean for each group
  group_means <- df %>%
    rowwise() %>%
    mutate(
      `EM-batch condition mix` = mean(c_across(all_of(groups$`EM-batch condition mix`)), na.rm = TRUE),
      `Denoised-batch mixing` = mean(c_across(all_of(groups$`Denoised-batch mixing`)), na.rm = TRUE),
      `Signal preservation` = mean(c_across(all_of(groups$`Signal preservation`)), na.rm = TRUE),
      `Data distortion` = mean(c_across(all_of(groups$`Data distortion`)), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(Method, `EM-batch condition mix`, `Denoised-batch mixing`, 
           `Signal preservation`, `Data distortion`)
  
  # Step 2: Apply 0-1 normalization (min-max scaling) to group means
  normalized_group_means <- group_means %>%
    mutate(across(
      -Method,  # Exclude the 'Method' column
      ~ (. - min(.)) / (max(.) - min(.)),  # Min-max normalization
      .names = "{col}_normalized"  # Add "_normalized" suffix
    ))
  
  # Step 3: Calculate total score for each method
  normalized_group_means <- normalized_group_means %>%
    rowwise() %>%
    mutate(
      Total_Score = mean(c_across(ends_with("_normalized")), na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Step 4: Prepare radar plot data using normalized columns
  radar_data <- normalized_group_means %>%
    select(ends_with("_normalized")) %>%
    as.data.frame()
  
  # Add max and min rows for scaling
  radar_data <- rbind(
    rep(1, ncol(radar_data)),  # Max values
    rep(0, ncol(radar_data)),  # Min values
    radar_data
  )
  
  # Assign column names for radar plot
  colnames(radar_data) <- gsub("_normalized", "", colnames(radar_data))  # Remove "_normalized" suffix
  
  # Step 5: Plot radar chart
  colors <- rainbow(nrow(group_means))  # Unique colors for each method
  
  par(mar = c(2, 2, 2, 2))  # Adjust margins
  print(radar_data)
  radar_data<-radar_data[,c(4,2,1,3)]
  # Rename the column
  colnames(radar_data)[colnames(radar_data) == "Data distortion"] <- "No data distortion"
  
  radarchart(
    radar_data,
    title=title,
    axistype = 1,
    pcol = colors,  # Line colors
    pfcol = alpha(colors, 0.4),  # Fill colors with transparency
    plwd = 2,  # Line width
    cglcol = "grey",  # Grid line color
    cglty = 1,  # Grid line type
    axislabcol = "black",  # Axis label color
    caxislabels = seq(0, 1, 0.2),  # Axis labels
    cglwd = 0.8,  # Grid line width
    vlcex = 1.5  # Increase vertex label size
  )
  
  # Add legend
  legend(
    "topright", legend = group_means$Method, col = colors, lty = 1, lwd = 2,
    cex = 1.2, bty = "n"  # Increase legend font size
  )
  
  # Step 6: Plot total scores for each method
  # Create a horizontal barplot with ggplot2
  ggplot(normalized_group_means, aes(x = reorder(Method, Total_Score), y = Total_Score, fill = Method)) +
    geom_bar(stat = "identity", width = 0.6) +  # Horizontal bars
    geom_text(aes(label = format(round(Total_Score, 2), nsmall = 2)),  # Add labels with 2 decimal places
              hjust = -0.2, size = 3.5, color = "black") +  # Adjust label position and size
    coord_flip() +  # Flip coordinates for horizontal bars
    labs(
      x = "Method",
      y = "Total Score"
    ) +
    scale_fill_manual(values = rainbow(nrow(normalized_group_means))[c(2,1,3)]) +  # Match colors to the radar plot
    theme_minimal() +  # Use a clean theme
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.position = "none"  # Remove legend
    )
}

df <- read.csv("./case_control/mouse.csv")
df<-df[,c(1,c(6:ncol(df)))]
generate_radar_plot(df,"Mouse")



