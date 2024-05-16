library(Seurat) 
library(ggplot2)
library(survival)
library(ggfortify)
library(AUCell)
library(dplyr)
library(readxl)
library(ggpubr)
library(tibble)
library(EnhancedVolcano)
library(S4Vectors)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(edgeR)
library(variancePartition)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL")


################################################
#Comparing relapse and initial samples 
################################################
regulons = read.table("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL/SCENICPlus/SCENICPlus_Pipeline/Output_BMP_TSpec_Final/eRegulon_filtered.tsv", header = T)
egrn_list = c()
for(egrn in unique(regulons$eRegulon_name)){
  egrn_list[[egrn]] = unique(regulons[regulons$eRegulon_name==egrn, ]$Gene)
}
regulon_sigs = egrn_list
names(egrn_list)

full_bulk = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/1_harmonizemtx_make_seurat/TPM.bulk.rna.umap.20.rds")

expr.matrix.full = full_bulk@assays$bulkRNA_TPM@counts
cells_rankings_full <- AUCell_buildRankings(expr.matrix.full, nCores = 12)
cells_AUC.Regulon.full <- AUCell_calcAUC(regulon_sigs, cells_rankings_full, 
                                         aucMaxRank=nrow(cells_rankings_full)*0.25, verbose = T)

full_bulk[["AUCell_Regulons"]] <- CreateAssayObject(cells_AUC.Regulon.full@assays@data$AUC)
rownames(full_bulk@assays$AUCell_Regulons)

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
bmp_sigs = list(BMP17, BMP_larger, T.spec_larger)
names(bmp_sigs) = c("BMP17", "BMPLarger", "TSpecLarger")

cells_AUC.BMP.full <- AUCell_calcAUC(bmp_sigs, cells_rankings_full, 
                                     aucMaxRank=nrow(cells_rankings_full)*0.25, verbose = T)
full_bulk[["AUCell_BMP"]] <- CreateAssayObject(cells_AUC.BMP.full@assays@data$AUC)
rownames(full_bulk@assays$AUCell_BMP)

full_bulk$USI %>% grep(pattern = "_R")
full_bulk$relapse.primary = "primary"
full_bulk$relapse.primary[grep(x = full_bulk$USI, pattern = "_R")] = "relapse"
table(full_bulk$relapse.primary)
full_bulk$patient_id = gsub("^(.*?)_.*$", "\\1", full_bulk$USI)
rownames(full_bulk@assays$AUCell_BMP)
patient_table = table(full_bulk$patient_id)
pairs = names(patient_table[patient_table>1])
pairs = c(pairs, paste0(pairs, "_R"))
length(pairs)

bmp.scores.with.mdata = cbind(full_bulk@assays$AUCell_BMP@counts %>% t() %>% as_tibble(), full_bulk@meta.data)
pairs %in% rownames(bmp.scores.with.mdata)
View(as.data.frame(bmp.scores.with.mdata$USI))
bmp.scores.with.mdata.paired <- bmp.scores.with.mdata[pairs, ]
table(bmp.scores.with.mdata.paired$relapse.primary)
table(bmp.scores.with.mdata.paired$patient_id)
table(bmp.scores.with.mdata.paired$ETP)


mean(na.omit(bmp.scores.with.mdata.paired$Percent.blasts.relapse.sample))
mean(na.omit(bmp.scores.with.mdata.paired$Percent.Blasts.Tumor.Sample.Diagnostic))

range(na.omit(bmp.scores.with.mdata.paired$Percent.blasts.relapse.sample))
range(na.omit(bmp.scores.with.mdata.paired$Percent.Blasts.Tumor.Sample.Diagnostic))


t.test(bmp.scores.with.mdata.paired$Percent.blasts.relapse.sample, 
       bmp.scores.with.mdata.paired$Percent.Blasts.Tumor.Sample.Diagnostic, paired = T)

bmp.scores.with.mdata.paired$percent.blasts.corrected = NULL

t.test(bmp.scores.with.mdata.paired$BMP17[1:27],
       bmp.scores.with.mdata.paired$BMP17[28:54], var.equal = F, paired = T)

pathways = rownames(full_bulk@assays$AUCell_BMP)
pathways <- gsub("_", "-", pathways)
for(pathway in pathways){
  print(pathway)
  p1 <- ggplot(bmp.scores.with.mdata.paired, 
               aes(x = relapse.primary, 
                   y = .data[[pathway]], 
                   fill = relapse.primary)) + 
    geom_boxplot() + geom_jitter(width = 0) + 
    theme_bw() +
    scale_fill_manual(values = c("purple", "darkgreen") %>% rev()) +
    theme(legend.position = "none") +
    labs(x = "", y = "") + geom_line(aes(group=patient_id)) +
    stat_compare_means(method = "t.test", paired = T)
  ggsave(p1, filename = paste0("Relapse_Figures/",pathway,".pdf"), width = 2.5, height = 4)
}

bmp.scores.with.mdata = cbind(full_bulk@assays$AUCell_Regulons@counts %>% t() %>% as_tibble(), full_bulk@meta.data)
bmp.scores.with.mdata.paired <- bmp.scores.with.mdata[pairs, ]
table(bmp.scores.with.mdata.paired$relapse.primary)
table(bmp.scores.with.mdata.paired$patient_id)
pathways = rownames(full_bulk@assays$AUCell_Regulons)
pathways <- gsub("_", "-", pathways)
for(pathway in pathways){
  print(pathway)
  p1 <- ggplot(bmp.scores.with.mdata.paired, 
               aes(x = relapse.primary, 
                   y = .data[[pathway]], 
                   fill = relapse.primary)) + 
    geom_boxplot() + geom_jitter(width = 0) + 
    theme_bw() +
    scale_fill_manual(values = c("purple", "darkgreen") %>% rev()) +
    theme(legend.position = "none") +
    labs(x = "", y = "") + geom_line(aes(group=patient_id)) +
    stat_compare_means(method = "t.test", paired = T)
  ggsave(p1, filename = paste0("Relapse_Figures/",strsplit(pathway, "/")[[1]][1],".pdf"), width = 2.5, height = 4)
}

##
#Adjusted for blast percentage
bmp.scores.with.mdata = cbind(full_bulk@assays$AUCell_BMP@counts %>% t() %>% as_tibble(), full_bulk@meta.data)
bmp.scores.with.mdata.paired <- bmp.scores.with.mdata[pairs, ]
bmp.scores.with.mdata.paired$percent.blasts.corrected[1:27] = bmp.scores.with.mdata.paired$Percent.Blasts.Tumor.Sample.Diagnostic[1:27]
bmp.scores.with.mdata.paired$percent.blasts.corrected[28:54] = bmp.scores.with.mdata.paired$Percent.blasts.relapse.sample[28:54]
bmp.scores.with.mdata.paired.with.blast.perc = subset(bmp.scores.with.mdata.paired, 
                                                      !(rownames(bmp.scores.with.mdata.paired) %in% c("PAUCIZ_R", "PAUFBJ_R", "PAUPTX_R", "PAVESF_R", 
                                                                                                      "PAUCIZ", "PAUFBJ", "PAUPTX", "PAVESF")))
pathways = rownames(full_bulk@assays$AUCell_BMP)
pathways <- gsub("_", "-", pathways)
for(pathway in pathways){
  print(pathway)
  p1 <- ggplot(bmp.scores.with.mdata.paired.with.blast.perc, 
               aes(x = relapse.primary, 
                   y = .data[[pathway]]/percent.blasts.corrected, 
                   fill = relapse.primary)) + 
    geom_boxplot() + geom_jitter(width = 0) + 
    theme_bw() +
    scale_fill_manual(values = c("purple", "darkgreen") %>% rev()) +
    theme(legend.position = "none") +
    labs(x = "", y = "") + geom_line(aes(group=patient_id)) +
    stat_compare_means(method = "t.test", paired = T)
  ggsave(p1, filename = paste0("Relapse_Figures/",pathway,"_Blast_Corrected.pdf"), width = 2.5, height = 4)
}

bmp.scores.with.mdata = cbind(full_bulk@assays$AUCell_Regulons@counts %>% t() %>% as_tibble(), full_bulk@meta.data)
bmp.scores.with.mdata.paired <- bmp.scores.with.mdata[pairs, ]
bmp.scores.with.mdata.paired$percent.blasts.corrected[1:27] = bmp.scores.with.mdata.paired$Percent.Blasts.Tumor.Sample.Diagnostic[1:27]
bmp.scores.with.mdata.paired$percent.blasts.corrected[28:54] = bmp.scores.with.mdata.paired$Percent.blasts.relapse.sample[28:54]
bmp.scores.with.mdata.paired.with.blast.perc = subset(bmp.scores.with.mdata.paired, 
                                                      !(rownames(bmp.scores.with.mdata.paired) %in% c("PAUCIZ_R", "PAUFBJ_R", "PAUPTX_R", "PAVESF_R", 
                                                                                                      "PAUCIZ", "PAUFBJ", "PAUPTX", "PAVESF")))
pathways = rownames(full_bulk@assays$AUCell_Regulons)
pathways <- gsub("_", "-", pathways)
for(pathway in pathways){
  print(pathway)
  p1 <- ggplot(bmp.scores.with.mdata.paired.with.blast.perc, 
               aes(x = relapse.primary, 
                   y = .data[[pathway]]/percent.blasts.corrected, 
                   fill = relapse.primary)) + 
    geom_boxplot() + geom_jitter(width = 0) + 
    theme_bw() +
    scale_fill_manual(values = c("purple", "darkgreen") %>% rev()) +
    theme(legend.position = "none") +
    labs(x = "", y = "") + geom_line(aes(group=patient_id)) +
    stat_compare_means(method = "wilcox.test", paired = T)
  ggsave(p1, filename = paste0("Relapse_Figures/",strsplit(pathway, "/")[[1]][1],"_Blast_Corrected.pdf"), width = 2.5, height = 4)
}



