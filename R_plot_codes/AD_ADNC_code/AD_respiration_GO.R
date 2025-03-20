library(clusterProfiler)
library(org.Hs.eg.db)

# Define the list of genes
gene_list <- c("ATP6", "COX2", "COX3", "CYTB", "GFAP", "GJA1", "LINC02742", "ND1", "ND2", "ND3", "ND4", "ND6", "PLCG2")

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check if there are unmapped genes
if (nrow(entrez_ids) < length(gene_list)) {
  cat("Some genes could not be mapped to Entrez IDs:\n")
  unmapped_genes <- setdiff(gene_list, entrez_ids$SYMBOL)
  print(unmapped_genes)
}

# Perform GO enrichment analysis
go_enrichment <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",       # Choose "BP", "CC", or "MF" for specific categories
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Display the results
head(go_enrichment)

# Visualize the results (e.g., barplot)
if (nrow(go_enrichment) > 0) {
  barplot(go_enrichment, showCategory = 10, title = "GO Enrichment Analysis")
} else {
  cat("No significant GO terms found.\n")
}