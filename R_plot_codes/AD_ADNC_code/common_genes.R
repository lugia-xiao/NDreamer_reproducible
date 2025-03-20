# Set environment language and working directory
Sys.setenv(LANG = "en")
setwd("D:/project/ndreamer_plot/")

# Read the dataframe
df <- read.csv("./AD_ADNC/down_significant_unparametric.csv", stringsAsFactors = FALSE)

# Filter to the top 100 rows (assuming rows represent the rank of significance)
df_top100 <- df[1:100, ]

# Get the number of cell types (columns)
num_cell_types <- ncol(df_top100)

# Find genes that appear in more than 80% of the cell types
gene_counts <- table(unlist(df_top100))  # Count occurrences of each gene across all columns
genes_80_percent <- names(gene_counts[gene_counts > 0.8 * num_cell_types])

# Output the genes
print("Genes appearing in more than 80% of the cell types among the top 100 most significant:")
print(genes_80_percent)
