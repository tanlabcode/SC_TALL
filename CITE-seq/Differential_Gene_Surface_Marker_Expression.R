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

####
#1
comparisons = list(c("T-specified-R", "T-specified-NR")) 

#2
comparisons = list(c("BMP-like-R", "BMP-like-NR")) 

#3
comparisons = list(c("T-specified-R", "BMP-like-NR")) 
#3b
comparisons = list(c("BMP-like-NR", "T-specified-R")) 

#4
comparisons = list(c("BMP-like-R", "T-specified-R")) 

#5
comparisons = list(c("BMP-like-R", "T-specified-NR")) 
####

group.by = c("comparison.bmp.vs.t.specified")
all.de.list= list()
human.tfs.and.cisbp = get.human.tfs()
tfs.intersect.10x = intersect(human.tfs.and.cisbp$human.tfs.and.cisbp, rownames(blasts.seurat@assays$RNA))

for (i in 1:length(comparisons)){
  
  group1 = comparisons[[i]][1]
  group2= comparisons[[i]][2]
  groupby = group.by[i]
  print(paste0(group1, "_vs_", group2))
  
  #tfs
  de.tfs = FindMarkers(blasts.seurat, 
                       verbose = T,
                       assay = "RNA",
                       features = tfs.intersect.10x,
                       logfc.threshold = 0, 
                       ident.1 = group1, 
                       ident.2 = group2, group.by = groupby, 
                       max.cells.per.ident = 1500)
  
  de.tfs$group = "de.tfs"
  de.tfs$gene = rownames(de.tfs)
  
  #genes
  
  de.genes = FindMarkers(blasts.seurat, 
                         verbose = T,
                         assay = "RNA", 
                         logfc.threshold = 0, 
                         ident.1 = group1, 
                         ident.2 = group2, 
                         group.by = groupby,
                         max.cells.per.ident = 1500)
  de.genes$group = "de.genes"
  de.genes$gene = rownames(de.genes)
  
  #cite markers
  de.surface = FindMarkers(blasts.seurat, 
                           verbose = T,
                           assay = "ADT", 
                           logfc.threshold = 0, 
                           ident.1 = group1, 
                           ident.2 = group2, 
                           group.by = groupby, 
                           max.cells.per.ident = 1500)
  de.surface$group = "de.surface"
  de.surface$gene = rownames(de.surface)
  
  #combine all
  all.de = rbind(de.tfs, de.genes, de.surface)
  
  all.de$comparison = paste0(group1, "_vs_", group2)
  all.de$enriched = ifelse(all.de$avg_log2FC > 0, yes = group1, no = group2)
  
  all.de.list[[paste0(group1, "_vs_", group2)]] = all.de
  filename2 = "all.de.list.rds"
  saveRDS(all.de.list, filename2)
  
}

#plot volcano plots
#filter surface markers

###
#1
de.all = all.de.list$`T-specified-R_vs_T-specified-NR`
de.all$comparison = "T-specified-R_vs_T-specified-NR"
de.all$enriched = ifelse(de.all$avg_log2FC < 0, yes = "T-specified-NR", no = "T-specified-R")

#2
de.all = all.de.list$`BMP-like-R_vs_BMP-like-NR`
de.all$comparison = "BMP-like-R_vs_BMP-like-NR"
de.all$enriched = ifelse(de.all$avg_log2FC < 0, yes = "BMP-like-NR", no = "BMP-like-R")

#3
de.all = all.de.list$`T-specified-R_vs_BMP-like-NR`
de.all$comparison = "T-specified-R_vs_BMP-like-NR"
de.all$enriched = ifelse(de.all$avg_log2FC < 0, yes = "BMP-like-NR", no = "T-specified-R")
#3b
de.all = all.de.list$`BMP-like-NR_vs_T-specified-R`
de.all$comparison = "BMP-like-NR_vs_T-specified-R"
de.all$enriched = ifelse(de.all$avg_log2FC < 0, yes = "T-specified-R", no = "BMP-like-NR")

#4
de.all = all.de.list$`BMP-like-R_vs_T-specified-R`
de.all$comparison = "BMP-like-R_vs_T-specified-R"
de.all$enriched = ifelse(de.all$avg_log2FC < 0, yes = "T-specified-R", no = "BMP-like-R")

#5
de.all = all.de.list$`BMP-like-R_vs_T-specified-NR`
de.all$comparison = "BMP-like-R_vs_T-specified-NR"
de.all$enriched = ifelse(de.all$avg_log2FC < 0, yes = "T-specified-NR", no = "BMP-like-R")
###


de.all = de.all %>% as_tibble() %>% filter(p_val_adj < 0.05)
surface.de = de.all %>% filter(group == "de.surface")
logcutoff = 1

de.plots= list()
de.plots$surface.loose = 
  ggplot(surface.de %>% filter(abs(avg_log2FC) > 0.25, p_val_adj < 1e-10), 
         aes(x = avg_log2FC, y = -log(p_val_adj, base = 20),
             color = enriched,
             label = gene)) + 
  geom_point(data = surface.de, color = "gray")+
  geom_point()+
  geom_text_repel(size = 5, max.overlaps = 10) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  theme(legend.position = "none")
de.plots$surface.loose

de.plots$surface.stringent = 
  ggplot(surface.de %>% filter(abs(avg_log2FC) > 0.5, p_val_adj < 1e-10), 
         aes(x = avg_log2FC, y = -log(p_val_adj, base = 20),
             color = enriched,
             label = gene)) + 
  geom_point(data = surface.de, color = "gray")+
  geom_point()+
  geom_text_repel(size = 5, max.overlaps = 10) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  theme(legend.position = "none")

de.plots$surface.stringent

gene.de = de.all %>% filter(group == "de.genes") %>% filter(!gene %in% human.tfs.and.cisbp$human.tfs.and.cisbp)

de.plots$gene.loose = ggplot(gene.de %>% filter(abs(avg_log2FC) > 0.4, p_val_adj < 1e-10), 
                             aes(x = avg_log2FC, y = -log(p_val_adj+1e-300, base = 10),
                                 color = enriched,
                                 label = gene)) + 
  geom_point(data = gene.de, color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 20) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  #scale_color_manual(values = c("#BE433D","#DDA164","#B89B74")) +
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  theme(legend.position = "none")
de.plots$gene.loose


de.plots$gene.stringent = ggplot(gene.de %>% filter(abs(avg_log2FC) > 0.5, p_val_adj < 1e-25), 
                                 aes(x = avg_log2FC, y = -log(p_val_adj+1e-300, base = 10),
                                     color = enriched,
                                     label = gene)) + 
  geom_point(data = gene.de, color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 20) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)") + #xlim(c(-3, 3))+ ylim(c(2, 300)) +
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  theme(legend.position = "none")
de.plots$gene.stringent

#make volcanos for tfs
tf.de = de.all %>% filter(group == "de.tfs") 
de.plots$tf.stringent = ggplot(tf.de %>% filter(abs(avg_log2FC) > 0.20, p_val_adj < 1e-05), 
                               aes(x = avg_log2FC, y = -log(p_val_adj+1e-300, base = 10),
                                   color = enriched,
                                   label = gene)) + 
  geom_point(data = tf.de, color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 20) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  theme(legend.position = "none")
de.plots$tf.stringent

de.plots$tf.loose = ggplot(tf.de %>% filter(abs(avg_log2FC) > 0.15, p_val_adj < 1e-03), 
                           aes(x = avg_log2FC, y = -log(p_val_adj+1e-300, base = 10),
                               color = enriched,
                               label = gene)) + 
  geom_point(data = tf.de, color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 20) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  theme(legend.position = "none")

de.plots$tf.loose

final.plot = cowplot::plot_grid(plotlist = de.plots, ncol = 3, byrow = F)


#Saving plots 
#1
pdf("Responder_Figures/TSpecified_R_vs_NR_de.plots.downsample1500.pdf", height = 4, width = 4)
de.plots
dev.off()

#2
pdf("Responder_Figures/BMP-like_R_vs_BMP-like-NR_de.plots.downsample1500.pdf", height = 4, width = 4)
de.plots
dev.off()

#3
pdf("Responder_Figures/TSpecified_R_vs_BMP-like-NR_de.plots.downsample1500.pdf", height = 4, width = 4)
de.plots
dev.off()
#3b
pdf("Responder_Figures/BMP-like-NR_vs_TSpecified_R_vs_de.plots.downsample1500.pdf", height = 4, width = 4)
de.plots
dev.off()

#4
pdf("Responder_Figures/BMP-like_R_vs_TSpecified_R_de.plots.downsample1500.pdf", height = 4, width = 4)
de.plots
dev.off()

#5
pdf("Responder_Figures/BMP-like_R_vs_TSpecified_NR_de.plots.downsample1500.pdf", height = 4, width = 4)
de.plots
dev.off()

