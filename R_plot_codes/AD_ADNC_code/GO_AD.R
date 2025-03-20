library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Set environment language and working directory
Sys.setenv(LANG = "en")
setwd("D:/project/ndreamer_plot/")

# Define the function for GO enrichment analysis and plotting
perform_go_enrichment <- function(df) {
  print(df)
  # Iterate over columns from the second to the last
  for (col_name in colnames(df)) {#colnames(df)[2:ncol(df)]
    print(col_name)
    # Extract gene list for the current column
    gene_list <- na.omit(df[[col_name]])
    print(gene_list)
    
    # Convert gene symbols to Entrez IDs
    gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    entrez_ids <- gene_entrez$ENTREZID
    
    # Perform GO enrichment analysis
    go_results <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL",  # Use ALL ontology
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    # Plot results
    if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
      print(paste("GO enrichment for:", col_name))
      print(dotplot(go_results, showCategory = 5) + ggtitle(paste("GO Enrichment for", col_name)))
    } else {
      print(paste("No significant GO terms found for:", col_name))
    }
  }
}

# Read the CSV file
df <- read.csv("./AD/up_significant_unparametric.csv")
df<-df[1:100,c("Lamp5","L2.3.IT")]#c("L2.3.IT","Oligodendrocyte","Microglia.PVM","Astrocyte","L4.IT")]

# Call the function
perform_go_enrichment(df)