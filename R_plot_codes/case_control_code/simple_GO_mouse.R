library(clusterProfiler)
library(org.Mm.eg.db)

# Define the list of genes
# ECCITE
#gene_list <- c('KLRD1', 'TERT', 'AC098973.2', 'FABP6', 'TFPI2', 'KLHL31', 'SCGB1B2P', 'NXPH3', 'CCDC160', 'LRCOL1', 'CD207', 'AL590708.2', 'CALML5', 'CD1C', 'NUTM2B', 'C21orf62', 'LINC00346', 'C22orf23', 'KRT19', 'BTNL2', 'RP11-159H22.2', 'SDK2', 'GUCY2F', 'ST18', 'RP4-777O23.1', 'HLA-DQA2', 'FAIM2', 'TNFRSF11B', 'KIRREL3', 'PGAM2', 'RP11-856M7.6', 'RP11-981G7.2', 'GDNF', 'CHRDL2', 'XCL1', 'KRTAP5-9', 'RP11-84D1.1', 'PLA2G4B', 'RP11-603J24.17', 'FAM203A', 'RP11-818F20.5', 'MGAT3', 'TAS2R46', 'OR13A1', 'GFRA3', 'DCN', 'PRRX1', 'AC115522.3', 'SCUBE3', 'FUT1', 'HTR2B', 'RND1', 'RP11-503N18.1', 'RP11-168G16.2', 'ZNF503-AS1', 'SLC10A1', 'FCGBP', 'MAPK4', 'RP11-651P23.5', 'EMILIN3', 'PDE6A', 'RP11-94L15.2', 'RP11-3D4.2', 'IFNL1', 'MGAT4C', 'LINC00284', 'VIPR2', 'C9orf163', 'RP11-700H6.1', 'TAL1', 'PCDHGB6', 'PLA2G2A', 'MSTN', 'RP5-901A4.1', 'AC004069.2', 'EBF4', 'EDN1', 'RP11-56B16.4', 'RP11-3D4.3', 'LOXL1', 'C4orf47', 'RP4-684O24.5', 'PDIA2', 'ACHE', 'CD300LG', 'C1orf173', 'ZNF667-AS1', 'RP11-113K21.4', 'RP11-55K13.1', 'TLL2', 'RP11-1020A11.2', 'RP11-394I13.3', 'RP11-546B15.1', 'SNCB', 'RP11-230C9.2', 'GULP1', 'APOBEC3A', 'PRR18', 'KCNQ2', 'CCDC147', 'RP11-192P3.4', 'ENTPD2', 'RP11-134K13.4', 'PLAC8L1', 'MMP1', 'AC109828.1', 'RP11-292E2.4', 'HYAL1', 'RP11-11C20.3', 'NTM', 'LINC00028', 'TSSK4', 'SNORD3B-1', 'SNAP25', 'RP11-727F15.12', 'RP3-439F8.1', 'CTD-2636A23.2', 'RP5-1098D14.1', 'LCN12', 'RP11-333I13.1', 'C17orf66', 'RP11-211C9.1', 'AL121578.2', 'RP5-902P8.12', 'CD300E', 'TTTY10', 'KRT23', 'CYR61', 'RP11-76E17.3', 'CNGA1', 'MTRNR2L2', 'SERPINB10', 'THEM5', 'LINC00449', 'SLC5A2', 'TACSTD2', 'CLEC4E', 'LAMC2', 'CTD-3065B20.2', 'RP4-799D16.1', 'NEXN', 'LINC00484', 'UNC5B-AS1', 'GPR17', 'RP6-109B7.5', 'RP11-158I3.3', 'IL1R2', 'RP11-524D16--A.3', 'ENPP2', 'AC007743.1', 'HR', 'RP11-20G13.3', 'RP11-74E22.3', 'EXTL1', 'CTD-2033A16.1', 'AC021218.2', 'RP11-20G13.2', 'COL22A1', 'HES4', 'RP11-244F12.3', 'EPHB1', 'HIST1H4I', 'CYP4F3', 'RP11-157D23.2', 'F2', 'RP3-455J7.4', 'FZD8', 'C9orf169', 'RP11-523H20.3', 'RP5-894D12.3', 'VGF', 'ZNF763.1', 'LXN', 'RP11-54O7.18', 'ODF3L1', 'CTC-575N7.1', 'GPR82', 'CHMP4C', 'RP11-317J10.2', 'MX2', 'RP11-481J2.2', 'BAAT', 'RASGRP1', 'BHLHE40-AS1', 'AC079466.1', 'DNAAF1', 'KLHL38', 'MCOLN2', 'ORM2', 'RP11-709B3.2', 'REM1', 'COL5A2', 'RP11-339B21.13', 'RP11-685M7.3', 'PCDHGA4', 'BACH2', 'RP11-122G18.8', 'RP11-389C8.3', 'USP18', 'RANBP3L')
# ASD1 
gene_list<-c('Ces2a', 'Enam', 'Gm15325', 'Gm12862', 'Gm29295', 'Gm6763', 'Krt4', 'B020011L13Rik', 'Noxa1', 'Usp17lc', '4930486F22Rik', 'Vmn1r21', 'Cd8a', 'Mat1a', 'Olfr1505', 'Gm26581', 'Gm15199', 'Doxl2', 'Nkx2-9', 'Gm7097', 'Tex21', 'Gm7298', '4930465M20Rik', 'Gm14135', 'Neurog3', 'Gm29017', 'Gh', 'Speer4d', '4930455H04Rik', 'Gm29554', 'Eppin', 'Gm26791', 'Tcf21', 'Tmem52b', '1700015E13Rik', 'Drd3', 'Gm42775', 'Gm11201', 'AY702103', 'Gm8334', 'Gm5416', 'Gm15945', 'Slc34a2', 'Ahsg', 'Gm13556', 'Clps', 'Gcnt3', 'Gpr87', '7630403G23Rik', 'Cd209a', 'Vmn2r10', 'Olfr1290', 'Zscan4c', 'Il24', 'Tex101', '2010106C02Rik', 'Gm37958', 'Gm12506', '1700109G15Rik', 'Gm8356', 'Zfp616', 'Gm20471', 'Entpd8', 'Pdc', '4930517E11Rik', 'Lpcat2b', 'Eddm3b', 'Gm13647', 'Pou4f3', 'Bcas3os2', 'Slco6d1', 'Cldn34b2', 'Gpr141', 'AF067063', 'Tmprss11e', 'Oas1f', 'Rln3', 'Gm4477', 'Them7', 'Gm13236', 'Cym', 'Try5', 'Cyp2ab1', 'Gm26816', 'Eif4e1b', 'Gm28178', 'Gm16081', 'Vmn2r82', 'Usp46os1', 'Gm15856', 'Lyg2', 'Fabp1', 'Gm14540', 'RP24-225H21.4', 'Oosp1', 'Gm1627', 'Gm26868', 'Gcsam', 'Gm15350', 'Gm5150', '4930533K18Rik', 'A930019D19Rik', 'Gm10722', 'Gm28111', '1700112H15Rik', '1700065J18Rik', 'AC132444.1', 'RP23-32A8.5', 'Apon', 'Ang5', 'Ceacam9', 'Clec12b', 'Igfn1', 'Gm18856', 'Olfr94', 'D830013O20Rik', 'Tespa1', 'Skint3', 'Otp', 'Hmx2', 'A530099J19Rik', 'Gm11730', 'Cd209d', '9830132P13Rik', 'Gm10804', 'Akr1d1', 'Crx', 'Olfr93', '4930599A14Rik', 'Esr2', 'Gm43400', 'Gm38404', 'Gm17197', 'Grk1', 'Tspear', 'RP23-278M8.1', 'Mettl7b', 'Gm9758', 'Gm15585', '4930583I09Rik', 'Vmn1r61', 'Spata31d1d', 'Cd160', 'Gm42439', 'Gm10339', 'Idi2', 'Gm27242', 'Gm5420', 'Gm42456', 'P2ry10', 'Nkg7', 'Pax8', '2700089I24Rik', 'S100g', '4933438A12Rik', 'Gm13580', 'Gm29246', 'Gm43253', '1700017N19Rik', 'Glp1r', 'Gimap3', 'Gstt4', 'Hist1h2br', 'Foxd3', 'Smco2', 'Dcpp1', '4930469G21Rik', 'Gm37446', 'Hoxa10', 'Gm27198', 'BC094916', 'Lcn3', 'Lhx8', 'Cpn2', 'Gm17173', '4930557J02Rik', 'Gm10718', 'Gm10324', '1700021F07Rik', 'Gm20753', 'Gm32250', '4930567H17Rik', 'Gpr132', 'Serpinb10', 'Usp44', 'Bpi', 'Gm831', 'Gm14372', 'Gm13008', 'Gm29392', 'Hoxd9', '1700023G09Rik', 'Scn10a', 'Gpr18', 'Gm11494', 'Gm26738', 'Slc51a', 'Crygs', 'Crxos', 'Gm26643')
dataset_name<-"ASD"

gene_list<-c('Ighv1-84', 'Ighv1-62-2', 'Ighv3-1', 'Igkv4-50', 'Tdgf1', 'Ppp2r2cos', '5730403I07Rik', 'Ighv11-2', '3632454L22Rik', 'Plxna4os2', 'Igkv4-74', 'Gm14471', 'Gm19466', 'Ighv1-42', 'Gm35147', 'Syndig1l', '1700095J12Rik', 'Igkv18-36', 'Gm32817', 'Ighv2-4', 'Igkv4-90', 'Gm9899', 'Gm26952', 'Ighv1-28', 'Vmn2r108', 'Gm19553', 'Slc22a20', 'Ighv9-2', 'Gm43409', '4930422M22Rik', 'Vmn1r201', 'Ighv7-3', 'Igkv8-21', 'Tlx2', 'Gm37027', 'Ly6h', 'Igkv4-80', 'Gm4598', 'Gm17122', 'Minar1', 'Gm36635', 'Trav22', 'Igkv8-19', 'Fhitos', 'Gm47922', 'Gkn1', 'Igkv4-68', 'Gm48942', 'Gm4473', 'Igkv11-125', 'Gm30146', 'Cd209f', 'Prps1l1', 'Vmn1r181', 'Fam163a', 'Scgb1a1', 'Gm50163', 'Igkv2-112', 'Ighv2-3', 'Ighv14-4', 'Csta3', 'Ighv1-61', 'Ighv1-4', 'Gm44850', 'Igkv4-62', 'Tmem262', 'Ighv3-8', 'Ighv1-5', 'Rsf1os1', 'Ttc21a', 'Gm13999', 'Gm42413', 'Prnd', 'Gm27007', 'Olfr675', 'Trav13-1', 'C330008A17Rik', 'Ighv1-7', 'Gm47047', 'Olfr1342', 'Chrna3', 'Gm11766', 'Gm11399', 'Gm32624', 'Adra2c', 'Akr1cl', 'Gm50362', 'Igkv12-98', 'Olfr398', 'Igkv3-12', 'Gm26646', 'Gm16070', 'Ighv1-34', 'Ighv14-3', 'Atp1b4', 'Tagln3', 'Olfm2', 'Gm48570', 'Gm48129', 'Gm34866', 'Gm49503', 'Fgf16', 'Comp', 'G530012D18Rik', 'Vmn2r110', 'Ctsg', 'Sprn', 'Cntn3', 'Ighv1-78', 'Igkv8-27', 'Dlx4os', 'Gm14798', 'Ighv9-3', 'Col28a1', 'Igkv7-33', 'Gm38403', 'Dsc3', '5330434G04Rik', 'Igkv4-63', 'Igkv4-79', 'Nmbr', 'Myt1l', 'Gm2800', 'Olfr9', 'Gm47350', 'Gm44696', 'Igkv6-32', 'Otogl', 'Igkv3-2', 'Gm31392', 'Ighv1-25', 'Apln', '4930444P10Rik', 'Kcna4', 'Kif6', 'Igkv13-84', 'Igkv2-116', 'Zbtb8b', 'Ighv5-15', 'Ighv5-16', 'Grm7', 'Ighv2-9', 'Trav8n-2', 'Gm42846', 'Gm44199', 'Gm42457', 'Gm12676', 'Rnf112', 'Iqub', 'Gm44238', 'Igkv6-17', 'Gm48677', 'Igkv5-48', 'Gm29554', 'Igkv1-133', 'Igkv4-77', 'Gm26676', 'Gm50222', '8030442B05Rik', 'Gm34408', 'Igkv6-25', '4930405A10Rik', 'Sfrp4', 'A430108G06Rik', 'Igll1', 'Ighv1-77', 'Cd209g', 'Gm32834', 'Ighv1-59', 'Gm13448', 'Igkv4-61', 'Iglv2', 'Lrrc15', 'Ighv2-9-1', 'Ighv1-20', 'Igkv3-5', 'Syt2', 'Klk10', 'Gm50255', 'Gm42982', 'Gm35438', 'Ctnna2', 'Them7', 'Gm3558', 'Ighv1-12', 'Gm47418', 'Folr2', 'Ighv5-17', 'Tac4', 'Krt4', 'Gm15056', 'Gm48843', 'Gcgr', 'Mgat5b', 'Gm4961', 'Draxin', 'Gm30108', '4930538E20Rik', 'Tuba3b', 'Igkv6-14')
dataset_name<-"Mouse"

print(length(gene_list))
gene_list<-gene_list[1:150]

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Check if there are unmapped genes
if (nrow(entrez_ids) < length(gene_list)) {
  cat("Some genes could not be mapped to Entrez IDs:\n")
  unmapped_genes <- setdiff(gene_list, entrez_ids$SYMBOL)
  print(unmapped_genes)
}

# Perform GO enrichment analysis
go_enrichment <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Mm.eg.db,
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
  barplot(go_enrichment, showCategory = 6, title = dataset_name)
} else {
  cat("No significant GO terms found.\n")
}