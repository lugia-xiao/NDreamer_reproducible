# Set environment language
Sys.setenv(LANG = "en")

# Set working directory
setwd("D:/project/ndreamer_plot/")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Define a function to create the plot
generate_plot <- function(input_df) {
  # Rename columns to correct issues
  colnames(input_df) <- colnames(input_df) %>% 
    gsub("X1.cLISI", "cLISI", .) %>% 
    gsub("kBET.Accept.Rate", "kBET", .) %>% 
    gsub("True Pos Rate", "True.Pos.Rate", .)
  print(input_df)
  # Define metric categories
  biological_conservation_metrics <- c("ASW_label", "ARI", "NMI", "cLISI")
  batch_correction_metrics <- c("bLISI_batch", "ASW_batch", "kBET_batch", "True.Pos.Rate_batch")
  batch_correction_perturbation_metrics <- c("bLISI_perturbation", "ASW_perturbation", 
                                             "kBET_perturbation", "True.Pos.Rate_perturbation")
  batch_correction_gene_target_metrics <- c("bLISI_gene_target", "ASW_gene_target", 
                                            "kBET_gene_target", "True.Pos.Rate_gene_target")
  #batch_correction_gene_target_metrics<-c()
  
  # Reshape data and calculate RelativeValue
  plot_data <- input_df %>%
    pivot_longer(-Method, names_to = "Metric", values_to = "Value") %>%
    filter(!is.na(Value), Metric %in% c(biological_conservation_metrics, 
                                        batch_correction_metrics, 
                                        batch_correction_perturbation_metrics, 
                                        batch_correction_gene_target_metrics)) %>%
    group_by(Metric) %>%
    mutate(RelativeValue = rescale(Value, to = c(0, 1))) %>%  # Scale values
    ungroup()
  
  # Calculate Total row for each method
  total_data <- plot_data %>%
    group_by(Method) %>%
    summarize(
      Metric = "Total",
      Value = mean(RelativeValue, na.rm = TRUE),
      RelativeValue = mean(RelativeValue, na.rm = TRUE)
    ) %>%
    mutate(Category = "Total")
  
  # Combine with the original data
  plot_data <- plot_data %>%
    mutate(Category = case_when(
      Metric %in% biological_conservation_metrics ~ "Biological Conservation",
      Metric %in% batch_correction_metrics ~ "Batch Correction",
      Metric %in% batch_correction_perturbation_metrics ~ "Condition Mixing - Perturbation",
      Metric %in% batch_correction_gene_target_metrics ~ "Condition Mixing - Gene Target"
    )) %>%
    bind_rows(total_data) %>%  # Add Total rows
    ungroup()
  
  # Reverse row order for plotting
  method_order <- plot_data %>%
    filter(Metric == "Total") %>%
    arrange(Value) %>% 
    pull(Method) %>% 
    rev()  # Reverse the order
  
  plot_data <- plot_data %>%
    mutate(Method = factor(Method, levels = rev(method_order)))
  print(plot_data)
  # Create the plot
  p <- ggplot(plot_data, aes(x = Metric, y = Method, fill = RelativeValue)) +
    geom_point(shape = 21, size = 19, color = "black") +  # Increase circle size
    geom_text(aes(label = round(Value, 3)), size = 5, vjust = 0.5, color = "#36395A") +  # Adjust text size and position, change text color
    scale_fill_gradientn(colors = c("#EFE7E5", "#CAD7EF")) +  # Custom color gradient
    scale_y_discrete(limits = levels(plot_data$Method)) +  # Respect factor levels
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14, color = "#36395A"),  # Bold, larger x-axis labels, custom color
      axis.text.y = element_text(face = "bold", size = 14, color = "#36395A"),  # Bold, larger y-axis labels (model names), custom color
      strip.text = element_text(size = 16, face = "bold", color = "#7D3432"),  # Increase size of category text, custom color for categories
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      fill = "Relative score"
    ) +
    facet_grid(~ Category, scales = "free_x", space = "free_x") +  # Separate by category
    geom_hline(yintercept = seq(1.5, n_distinct(plot_data$Method) - 0.5, by = 1), 
               linetype = "dashed", color = "grey")  # Add dashed lines between rows
  
  print(p)
}


# Read the data
paths<-c("./experimental/ASD1.csv","./experimental/ECCITE.csv")

for (pathi in paths){
  generate_plot(read.csv(pathi))
}