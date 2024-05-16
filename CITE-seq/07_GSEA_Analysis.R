source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/JX_scFunctions.R')

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL")
###############################
#Response Non Responder
##############################
blasts.seurat = readRDS("SCENICPlus/etp.25.g1.downsample.1711.each.patient.rds")
blasts.seurat.ann.metadata = blasts.seurat@meta.data %>% 
  mutate(comparison.bmp.vs.t.specified = case_when(response.binary == "mrd.over.01" & predicted.cell.type.short %in% c("HSPC", "LMPP", "CLP", "ETP") & ETP == "ETP" ~ "BMP-like-NR",
                                                   response.binary == "no.mrd" & predicted.cell.type.short %in% c("Pro-T", "Pre-T") & ETP == "ETP" ~ "T-specified-R", 
                                                   response.binary == "no.mrd" & predicted.cell.type.short %in% c("HSPC", "LMPP", "CLP", "ETP") & ETP == "ETP" ~ "BMP-like-R",
                                                   response.binary == "mrd.over.01" & predicted.cell.type.short %in% c("Pro-T", "Pre-T") & ETP == "ETP" ~ "T-specified-NR"))

blasts.seurat = AddMetaData(blasts.seurat, metadata = blasts.seurat.ann.metadata)

table(blasts.seurat$comparison.bmp.vs.t.specified)
sum(table(blasts.seurat$comparison.bmp.vs.t.specified))
dim(blasts.seurat)

Idents(blasts.seurat) = "comparison.bmp.vs.t.specified"

comparisons = list(c("T-specified-R", "T-specified-NR"),
                   c("BMP-like-R", "BMP-like-NR"),
                   c("BMP-like-NR", "T-specified-R"),
                   c("BMP-like-R", "T-specified-R"),
                   c("BMP-like-R", "T-specified-NR"),
                   c("BMP-like-NR", "T-specified-NR") )
groupby = "comparison.bmp.vs.t.specified"

#Run FindMarkers with Seurat v4 or v5
for(diff in comparisons){
  print(diff)
  group1 = diff[1]
  group2 = diff[2]
  de.genes = FindMarkers(blasts.seurat, min.pct = 0.001, 
                         verbose = T,
                         assay = "RNA", 
                         logfc.threshold = 0, 
                         ident.1 = group1, 
                         ident.2 = group2, 
                         group.by = groupby,
                         max.cells.per.ident = 1500)
  de.genes$group = "de.genes"
  de.genes$gene = rownames(de.genes)
  saveRDS(de.genes, paste0("Responder_Figures/",group1,"_vs_",group2,"_All_Genes_v5.rds"))
}


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

#Prepare GSEA
library(S4Vectors)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(gridExtra)

pathways.hallmark <- msigdbr(species = "human", category = "H")
pathway_list = split(x = pathways.hallmark$gene_symbol, f = pathways.hallmark$gs_name)
pathway_list[["BMP17"]] <- BMP17
pathway_list[["BMP_larger"]] <- BMP_larger
pathway_list[["T.spec_larger"]] <- T.spec_larger

#Run GSEA
for(diff in comparisons){
  print(diff)
  group1 = diff[1]
  group2 = diff[2]
  de.genes = readRDS(paste0("Responder_Figures/",group1,"_vs_",group2,"_All_Genes_v5.rds"))
  res <- de.genes %>% dplyr::select(gene, avg_log2FC) 
  ranks <- deframe(res)
  
  GSEA_result <-  fgsea(pathways = pathway_list, stats = ranks)
  write.table(as.data.frame(GSEA_result[,-"leadingEdge"]), file = paste0("Responder_Figures/",group1,"_vs_",group2,"_GSEA_Barplot.txt"), sep = '\t', quote = F)
  
  sig_pathways <- GSEA_result[GSEA_result$padj<0.05, ]
  sig_pathways <- sig_pathways[order(sig_pathways$NES), ]
  sig_pathways$pathway <- factor(
    sig_pathways$pathway, 
    levels = sig_pathways$pathway[order(sig_pathways$NES)]
  )
  top_pathways <- head(sig_pathways[sig_pathways$NES<0, ], 10)
  bottom_pathways <- tail(sig_pathways[sig_pathways$NES>0, ], 10)
  sig_pathways_filter <- rbind(top_pathways, bottom_pathways)
  sig_pathways_filter$pathway <- factor(
    sig_pathways_filter$pathway,
    levels = sig_pathways_filter$pathway[order(sig_pathways_filter$NES)]
  )
  p1 <- ggplot(sig_pathways_filter, aes(x = pathway, y = NES, fill = -log10(padj))) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 7.5) +
    labs(x = "", y = "Normalized Enrichment Score", fill = "FDR") +  # Update axis and legend labels
    theme_bw() + coord_flip()+
    theme(
      axis.text.x = element_text(vjust = 0, hjust = 0.5, family = "sans", size = 8, color = "black"),
      axis.text.y = element_text(family = "sans", size = 8, color = "black"),
      axis.line = element_blank(),
      plot.margin = margin(l = 5, r = 5, b = 5, t = 5, unit = "pt"), 
      panel.grid = element_blank())
  p1
  ggsave(paste0("Responder_Figures/",group1,"_vs_",group2,"_GSEA_Barplot.pdf"), 
         plot = p1, device = 'pdf', height = 2, width = 6)
  
  p3a <- plotEnrichment(pathway_list[["BMP17"]],
                        ranks) + labs(title="BMP17") + theme_bw()
  p3b <- plotEnrichment(pathway_list[["BMP_larger"]],
                        ranks) + labs(title="BMP_larger") + theme_bw()
  p3c <- plotEnrichment(pathway_list[["T.spec_larger"]],
                        ranks) + labs(title="T.spec_larger") + theme_bw()
  p6 <- grid.arrange(p3a, p3b, p3c, nrow = 1)
  p6
  ggsave(paste0("Responder_Figures/",group1,"_vs_",group2,"_GSEA_3Enrich.pdf"), 
         plot = p6, device = 'pdf', height = 1.7, width = 6.7)         
  
}
the p[at]
