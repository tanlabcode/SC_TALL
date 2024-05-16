library(Seurat)
library(survival)
library(dplyr)
library(limma)
library(edgeR)
library(ggpubr)
library(gridExtra)
library(DGEobj.utils)
library(AUCell)
library(Matrix)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/T_ALL/BulkRNAseq")

load("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/data_sharing/TALL_X01_normalized_1362.Rdata")
lookup_table <- setNames(gene_anno$geneSymbol, gene_anno$geneID)
rownames(filtered.counts.corrected) <- lookup_table[rownames(filtered.counts.corrected)]
featureLength = gene_anno[gene_anno$geneSymbol %in% rownames(filtered.counts.corrected), ]
featureLength <- featureLength[!duplicated(featureLength$geneSymbol), ]
filtered.counts.corrected = filtered.counts.corrected[featureLength$geneSymbol, ]

dim(filtered.counts.corrected)
dim(featureLength)
identical(rownames(filtered.counts.corrected), featureLength$geneSymbol)
lengths = featureLength$length

fpkm_new = convertCounts(filtered.counts.corrected, unit="FPKM", geneLength=lengths)

#Previous drug data 
fpkm = as.data.frame(read.csv("pharmacotyping_ped_rnaseq_fpkm_ALLids_0823.csv"))
rownames(fpkm) = fpkm$GeneName
fpkm$GeneID = NULL
fpkm$GeneName = NULL
rownames(fpkm)[rownames(fpkm) == "KIAA0125"] <- "FAM30A"

common_genes = intersect(rownames(fpkm), rownames(fpkm_new))
fpkm_common = fpkm[common_genes, ]
fpkm_new_common = fpkm_new[common_genes, ]

fpkm_all = merge(fpkm_new_common, fpkm_common, by = "row.names", all = TRUE)
dim(fpkm_common)
dim(fpkm_new_common)
dim(fpkm_all)
rownames(fpkm_all) = fpkm_all$Row.names
fpkm_all$Row.names = NULL

#BMP17 Z-scores
BMP17 = c("S100A4", "LGALS1", "FAM30A", "IGLL1", "MEF2C", 
          "CTSW", "PRDX1", "HCST", "HSH2D", "KLF2", "VAMP8", 
          "HOPX", "CYBA", "MT-ND3", "ENO1", "PEBP1", "CD44")

BMP_larger <- c("PRDX1", "S100A4", "BCL11A", "RAMP1", "SPINK2", "IGFBP7", "CD44", "C1QTNF4", "DNAJC1", "COMMD3", 
                "PRSS57", "IGLL1", "LGALS1", "GYPC", "HCST", "ITGA4", "CTSW", "KIAA0087", "CYBA", "MEF2C", 
                "MATK", "IFITM2", "CYTL1", "BAALC", "SMIM24", "HOXA10", "PTH2", "TENT5A", "CAPG", "PEBP1", 
                "HSH2D", "VAMP8", "LPCAT2", "EEF1A2", "HOXA9", "FAM30A", "BMI1", "GOLGA8N", "DDAH2", "MT-ND3", 
                "MT-ND4L", "MCTP2", "ENO1", "AHNAK", "TRH", "XIST", "CEBPG", "IGKC", "HOPX", "GZMA", "MSI2", 
                "GSTP1", "MTRNR2L8", "APLP2", "PTPRE", "ZNF503", "AC103591.3", "KLF2", "DAD1", "CST3", "PLCG2", 
                "TRPS1", "PRDX2", "ARL4C", "MDK", "ACY3")
T.spec_larger  <- c("HES4", "MAL", "SCGB3A1", "BCL11B", "BTG3", "TASP1", "CHI3L2", "TRGC2", "PMEPA1", "IFITM3",
                    "NDFIP1", "SDCBP", "BEX3", "TCF7", "SH3BP5", "LAT", "MYO7B", "ITGA1", "CD3E", "HERPUD1",
                    "CRYBG1", "PCGF5", "PMAIP1", "CD99", "MIR4435-2HG", "MZB1", "SIT1", "LY6H", "MGAT4A",
                    "PRSS2", "L1TD1", "LCK", "ITM2B", "UBE2B", "CLDN5", "CCND3", "LPAR6", "SNTG2", "YBX3",
                    "PCDH10", "LTB", "MALT1", "TSPYL2", "SAT1", "DNTT", "RBM38", "PPP1CB", "MAP1A", "SDF2L1",
                    "DDIT4", "CYTOR", "TUBB2A", "HSP90B1")

fpkm_log2 = log2(fpkm_all+0.1)
dim(fpkm_log2)

drugs.to.keep = c("Venetoclax", "Ibrutinib", "Prednisolone", "Navitoclax", 
                  "Nelarabine", "Mercaptopurine", "Vincristine", 
                  "Daunorubicin", "Clofarabine")

# drugs.to.keep = c("Asparaginase", "Bortezomib", "CHZ868", "Cytarabine", "Dasatinib", "Dexamethasone", "Panobinostat", 
#                   "Ruxolitinib", "Thioguanine", "Trametinib", "Vorinostat")

signatureScore_zscore_robust <- function(matrix, pfeatures, nfeatures = NULL,
                                         vars.to.regress = NULL,
                                         score.name = 'Apop'){
  pfeatures = pfeatures[pfeatures %in% rownames(matrix)]
  mtx = matrix
  rzscore <- function(x){
    sd.robust = sd(x)/2 + mad(x)/2
    x.robust = (x - median(x))/sd.robust
    if(sd.robust <= 0) x.robust = x - median(x)
    return(x.robust)
  }
  pscore = mtx[pfeatures, ,drop = F]
  pscore = t(apply(pscore, 1, rzscore))
  pscore = Matrix::colSums(pscore)
  if(length(nfeatures) > 0) {
    nscore = mtx[nfeatures, ,drop = F]
    nscore = na.omit(nscore)
    nscore = t(apply(nscore, 1, rzscore))
    nscore = Matrix::colSums(nscore)
    scores = (pscore - nscore)/length(c(pfeatures, nfeatures))
  }else{
    scores = pscore/length(pfeatures)
  }
  return(scores)
}

z_scores = signatureScore_zscore_robust(fpkm_log2, pfeatures = BMP_larger, nfeatures = T.spec_larger, score.name = "BMP119")

kmeans_result <- kmeans(z_scores, centers = 2, nstart = 100, iter.max=1000)
cluster_means <- tapply(z_scores, kmeans_result$cluster, mean)
ordered_levels <- order(cluster_means)
levels <- factor(kmeans_result$cluster, levels = ordered_levels)
table(levels)

z_levels = levels
names(z_levels) = names(z_scores)

plots = c()
plots2 = c()
for(name in drugs.to.keep){
  drug_mdata = read.table("Drug_Concentrations_All.tsv", sep = '\t', header = T)

  drug_mdata = drug_mdata[drug_mdata$Molecular.subtype %in% c("T-ALL", "ETP"), ]
  drug_mdata = drug_mdata[!is.na(drug_mdata[[name]]), ]
  
  if(nrow(drug_mdata)>2){
    rownames(drug_mdata) = drug_mdata$Patient.ID
    drug_cases = drug_mdata$Patient.ID
    dim(drug_mdata)
    
    print(name)
    
    existing_columns <- intersect(drug_cases, colnames(fpkm_all))
    
    drug_mdata_filter = drug_mdata[existing_columns, ]
    drug_mdata_filter$Dataset = factor(drug_mdata_filter$Dataset)
    drug_mdata_filter$BMP119_Score <- as.numeric(z_scores[existing_columns])
    drug_mdata_filter$BMP119_Level <- factor(as.numeric(z_levels[existing_columns]))
    
    # Create  plots
    drug_mdata_filter$drug_name = drug_mdata_filter[[name]]
    correlation <- cor.test(-log2(drug_mdata_filter$drug_name), drug_mdata_filter$BMP119_Score, method = "spearman")
    print(correlation)
    print(dim(drug_mdata))
    p1 <- ggplot(drug_mdata_filter, aes(x = BMP119_Score, y = -log2(drug_name), color = Molecular.subtype, shape = Dataset)) +
      geom_point() +  # Add points for each data point
      geom_smooth(aes(group = 1), method = "lm", se = F, color = "blue") +  # Add linear regression line
      theme_bw() +  # Use a white and gray theme
      labs(y = paste0("-log2(",name," Sensitivity)"), x = "BMP119_Score", title = name) +  
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      guides(color = FALSE, shape = F) 
    cor_text <- paste("rho =", round(correlation$estimate, 2), "\np =", signif(correlation$p.value, digits = 2))
    p1 <- p1 +
      annotate("text", y = min(na.omit(-log2(drug_mdata_filter$drug_name))), x = min(drug_mdata_filter$BMP119_Score),
               label = cor_text, hjust = -1.1, vjust = -0.5, size = 4, color = "black")
    
    plots[[name]] = p1
    
    p2 = ggplot(drug_mdata_filter, aes(x = BMP119_Level, y = log2(drug_name), fill = BMP119_Level)) +
      geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + stat_compare_means(method = "wilcox.test") + 
      labs(x = "BMP119 Level", y = "log2(Drug Sensitivity)", title = name) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      guides(fill = FALSE) +
      scale_x_discrete(labels = c("Low", "High"))
    plots2[[name]] = p2
  }
  
}
plots = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                     plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], nrow = 2)
plots2 = grid.arrange(plots2[[1]], plots2[[2]], plots2[[3]], plots2[[4]],
                      plots2[[5]], plots2[[6]], plots2[[7]], plots2[[8]], plots2[[9]], nrow = 2)
ggsave(plots, file = "Figures/Selected_Correlations.pdf", width = 12, height = 5)
ggsave(plots2, file = "Figures/Selected_Correlations_Boxplot.pdf", width = 10, height = 8)

# plots = grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
#                      plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]],
#                      plots[[10]], plots[[11]], nrow = 3)
# ggsave(plots, file = "Figures/Other_Selected_Correlations.pdf", width = 12, height = 5)
# 



