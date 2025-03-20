gene_list1<-c('FTH1', 'PCDH9-AS2', 'FTL', 'P2RY14', 'IPO9-AS1', 'ANGPT2', 'FN1', 'SHOC1', 'PDE8A', 'HHIP', 'LINC02899', 'TIMP3', 'B3GNT5', 'ENO4', 'PPP1R9A-AS1', 'SLC7A14-AS1', 'PDE7B-AS1', 'TMEM165', 'RFX4', 'CFTR', 'CTNNA3', 'BACH1', 'LINC02930', 'ABCB1', 'FGF1', 'CPB2', 'IL18', 'TGFB1', 'MARCKSL1', 'CCK', 'HIF1A-AS3', 'CRH', 'LPAR6', 'TP53TG5', 'DLEU1', 'EYS', 'SST', 'STARD13-AS', 'PACRG-AS3', 'TAC1', 'GNB4', 'GNG12-AS1', 'DBNDD2', 'MIR325HG', 'LOC105369147', 'TTN', 'HLA-A', 'LOC613266', 'ZBTB20-AS5', 'DPP10-AS3', 'MKRN3', 'CXCL14', 'KCNQ1', 'PRKCH', 'NEDD9', 'A2M', 'TRPC3', 'SAT1', 'VIP', 'CBR1', 'TFRC', 'ZNF385D-AS2', 'RBPMS', 'LOC105378311', 'SHISA8', 'CMTM8', 'CST3', 'RIN3', 'DISC1FP1', 'FAM149A', 'PKP2', 'SLC16A1', 'GAD2', 'MT3', 'CP', 'PCDH11Y', 'SLAIN1', 'RUBCNL', 'ATP1B2', 'WIF1', 'NXPH2', 'KIF26B', 'SLC25A13', 'LINC02694', 'EMC10', 'KMO', 'RNASET2', 'LHFPL3-AS1', 'FAM222A-AS1', 'LINC01060', 'HAPLN2', 'BCAN', 'ANKRD37', 'MGST1', 'LMCD1-AS1', 'LOC105371366', 'TAC3', 'ADAM12', 'SVEP1', 'SLC26A5-AS1')
gene_list2<-c('PCDH9-AS2', 'P2RY14', 'FTH1', 'PPP1R9A-AS1', 'TIMP3', 'HHIP', 'DPYD-AS1', 'CTNNA3', 'CHST9', 'ANGPT2', 'B3GNT5', 'ABCB1', 'SHOC1', 'RFX4', 'IPO9-AS1', 'DLEU1', 'FN1', 'DENND3', 'P2RY12', 'PCED1B', 'SLC7A14-AS1', 'EYS', 'LOC105369147', 'PDE7B-AS1', 'TP53TG5', 'TMEM165', 'MGST1', 'LY86', 'LINC01877', 'LINC02930', 'TGFB1', 'PDE4B', 'LINC02694', 'LPAR6', 'UST', 'PALM2AKAP2', 'CFTR', 'IL18', 'RIN3', 'HIF1A-AS3', 'PRKCH', 'ENO4', 'CFAP418-AS1', 'LINC02899', 'RBPMS', 'PACRG-AS3', 'RUNX2', 'LOC613266', 'SHISA6', 'NEDD9', 'FGF1', 'STARD13-AS', 'CPB2', 'FTL', 'SAT1', 'TTN', 'FAM149A', 'SRGAP2B', 'BMP2K', 'KMO', 'ABCC4', 'FPR3', 'ZNF804B', 'LINC00320', 'A2M', 'CCK', 'ACSS3', 'SRGAP2', 'LMCD1-AS1', 'LINC01090', 'SELE', 'MYO1F', 'DPP10-AS3', 'LOC105370315', 'GASK1B-AS1', 'PDE8A', 'TRAF3IP3', 'DBNDD2', 'MAN2A1', 'SCHLAP1', 'RASGRP3', 'NPAS3', 'TFRC', 'CP', 'ADAM12', 'LDB3', 'SLC16A1', 'HAPLN2', 'LOC100506474', 'NTN1', 'CIITA', 'PAMR1', 'CPQ', 'PCDH11Y', 'CHRM5', 'LOC107985126', 'TACR1', 'CEBPD', 'LYVE1', 'CBR1')

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Set environment language and working directory
Sys.setenv(LANG = "en")
setwd("D:/project/ndreamer_plot/")

# Define the function for GO enrichment analysis and plotting
perform_go_enrichment <- function(df) {
  # Iterate over columns from the second to the last
  for (col_name in colnames(df)[1:ncol(df)]) {#colnames(df)[2:ncol(df)]
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
      print(dotplot(go_results, showCategory = 9) + ggtitle(paste("GO Enrichment for", col_name)))
    } else {
      print(paste("No significant GO terms found for:", col_name))
    }
  }
}

# Read the CSV file
df <- read.csv("./AD/down_significant_parametric.csv")
df['GABAergic']<-gene_list1
df['Glutamatergic']<-gene_list2
df<-df[,c('GABAergic','Glutamatergic')]

# Call the function
perform_go_enrichment(df)