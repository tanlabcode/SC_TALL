library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL")

rna.all = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/60_v1_figures_ETPDataset/Objects/t.all.40.with.cMeta.rds")

DefaultAssay(rna.all) = "RNA"
table(rna.all$is.blast.viscello)
table(rna.all$sample.group.short)
DimPlot(rna.all, group.by = "level.1.anno.viscello", reduction = "umap_50_defaultVEG") + coord_fixed()
length(table(rna.all$orig.ident))

Idents(rna.all) = "orig.ident"
rna.all[["RNA"]] <- split(rna.all[["RNA"]], f = rna.all$orig.ident)
rna.all <- NormalizeData(rna.all)
rna.all <- FindVariableFeatures(rna.all)
rna.all <- ScaleData(rna.all)
rna.all <- RunPCA(rna.all)
rna.all <- RunUMAP(rna.all, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
saveRDS(rna.all, "TALL_40_Unintegrated_JS.rds")

rna.all <- IntegrateLayers(
  object = rna.all, method = RPCAIntegration,  
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE
)
rna.all <- RunUMAP(rna.all, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
#saveRDS(rna.all, "TALL_40_RPCA_Integrated_JS.rds")

colors = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure2_R_vs_NR/_colors_cell_type_short.rna-atac.rds")

rna.all$level.1.anno.viscello.sample.group.short = factor(rna.all$level.1.anno.viscello.sample.group.short, 
                                                          levels = c("Healthy_B", "Healthy_Myeloid", "Healthy_Progenitor", 
                                                                     "Healthy_T/NK", "ETP_Blast", "Near-ETP_Blast", "Non-ETP_Blast"))
annotation_colors = c("#A6CEE3", "#99CD91", "#A264DA", "#459DA5", "#BE433D", "#DDA164", "#B89B74")
p1 = DimPlot(rna.all, group.by = "level.1.anno.viscello.sample.group.short", reduction = "umap.rpca",
             raster = T, raster.dpi = c(3000,3000), pt.size = 3, alpha = 0.7, shuffle = T, cols = annotation_colors) + coord_fixed()
ggsave(p1, file = "Annotation_Figures/Integrated_UMAP_Annotation.pdf", width = 6, height = 6)

p2 = DimPlot(rna.all, group.by = "orig.ident", reduction = "umap.rpca",
             raster = T, raster.dpi = c(3000,3000), pt.size = 3, alpha = 0.7, shuffle = T) + coord_fixed()
ggsave(p2, file = "Annotation_Figures/Integrated_UMAP_Sample.pdf", width = 9, height = 9)

#Correlation with pathology
load("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/data_sharing/TALL_X01_normalized_1362.Rdata")

rna.all = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/60_v1_figures_ETPDataset/Objects/t.all.40.with.cMeta.rds")
atac.all = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/21_label_transfer_from_rna_ref/atac.ref.after.label.transfer.harmonized.labels.fix.b.mislabeled.with.maturity.rds")

bulk.rna.object.umap20 <- readRDS("/mnt/isilon/tan_lab/sussmanj/Temp/T_ALL/BulkRNAseq/vst.no.relapse.bulk.rna.umap.20.rds")

bulk_blast_prop = data.frame(Bulk_Blast = bulk.rna.object.umap20$Percent.Blasts.Tumor.Sample.Diagnostic, Sample.Type = bulk.rna.object.umap20$Tumor.Specimen.Type)
rna_blast_prop = prop.table(table(rna.all$orig.ident, rna.all$is.blast.viscello), margin = 1)
colnames(rna_blast_prop) = c("RNA_Healthy", "RNA_Blast")
patients = unique(rna.all$orig.ident)

rna_blast_prop = rna_blast_prop[patients, ]
bulk_blast_prop = bulk_blast_prop[patients, ]

#rna_blast_prop = prop.table(table(atac.all$orig.ident, atac.all$Ctype_binned), margin = 1)
#names(rna_blast_prop) = c("RNA_Healthy", "RNA_Blast")

merged_df <- as.data.frame.matrix(rna_blast_prop)
merged_df$Bulk_Blast = as.numeric(bulk_blast_prop$Bulk_Blast)
merged_df$Sample.Type = bulk_blast_prop$Sample.Type
merged_df = na.omit(merged_df)
merged_df$Bulk_RNA_Diff = merged_df$Bulk_Blast - merged_df$RNA_Blast
merged_df$Abs_Bulk_RNA_Diff = abs(merged_df$Bulk_RNA_Diff)
merged_df$All = "All"

#merged_df = merged_df[merged_df$Sample.Type == "bone marrow", ]

cor.test(merged_df$Bulk_Blast, merged_df$RNA_Blast)
ggplot(merged_df, aes(x = Bulk_Blast, y = RNA_Blast)) +
  geom_point(color = "black") +  # Scatter points
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Linear regression line
  labs(
    x = "Pathology Estimate",
    y = "scRNA-seq"
  ) + theme_bw() 

mean(merged_df$Bulk_RNA_Diff)
ggplot(merged_df, aes(x = All, y = Bulk_RNA_Diff)) +
  geom_boxplot() +   
  geom_jitter(width = 0.2) + theme_bw() + 
  labs(y = "Bulk RNA Diff") 

mean(merged_df$Abs_Bulk_RNA_Diff)
ggplot(merged_df, aes(x = All, y = Abs_Bulk_RNA_Diff)) +
  geom_boxplot() +   
  geom_jitter(width = 0.2) + theme_bw() + 
  labs(y = "Bulk RNA Diff") 

hist(c(merged_df$RNA_Blast, merged_df$Bulk_Blast))
merged_df$patient = rownames(merged_df)
ggpaired(merged_df, cond1 = "Bulk_Blast", cond2 = "RNA_Blast", id = "patient") + 
  stat_compare_means(paired = T)
wilcox.test(merged_df$Bulk_Blast, merged_df$RNA_Blast, paired = T)
dim(merged_df)
