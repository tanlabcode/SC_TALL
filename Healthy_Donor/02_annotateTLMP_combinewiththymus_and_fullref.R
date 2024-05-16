

thymus = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/ReferenceMapping/proj_Thymus6Combined/thymus6combined_with_proj_anno.RDS")

full.ref = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/other_r_objects/tan_lab_combined_ref/4_23_full_ref_DN_early_late.RDS")

tlmps = readRDS("/mnt/isilon/tan_lab/xuj5/HealthyReference/SeuratObjects/Thymus/TLMP_4M_RNASeq.rds")

FeaturePlot(tlmps, "IL7R")
DimPlot(tlmps ,label = T)
tlmps@meta.data = tlmps@meta.data %>% rename(RNA_snn_res.1 = "TLMP_RNA_snn_res.1")


#do degs for all clusters
markers = FindMarkers(tlmps, ident.1 = 8, ident.2 = 5,slot = "counts", logfc.threshold = 0.5, latent.vars = "nCount_RNA")

batch = markers %>% as_tibble(rownames = "gene") %>% filter(avg_logFC > 0)


#setup function to runumap
run.umap.10.20 = function(object){
  
  hello = object %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress = c('perc.mito', 'nCount_RNA', 'S.Score', 'G2M.Score')) %>% 
    RunPCA(npcs = 30) %>% RunUMAP(dims = 1:20, reduction.name = "umap20") %>% RunUMAP(dims = 1:10, reduction.name = "umap10") %>% FindNeighbors() %>% FindClusters(res = 0.3)
  
  return(hello)
}

#dn.only
dn.only = subset(full.ref, combined.celltype %in% c("DN(Q)", "DN(P)", "DN(early)"))
all.dns = merge(dn.only, tlmps) %>% run.umap.10.20()


DimPlot(all.dns, reduction = "umap20", label = F, group.by = "sample.name")
DimPlot(all.dns, reduction = "umap20", label = T, group.by = "combined.celltype")
DimPlot(all.dns, reduction = "umap20", label = T, group.by = "TLMP_RNA_snn_res.1")
FeaturePlot(all.dns, c("IGHM", "CD74"))
all.dns@meta.data$RNA_snn_res.0.3
all.dns@meta.data$cell.name = rownames(all.dns@meta.data)

FeaturePlot(all.dns, "CD1A", reduction = "umap20")
big.ref = merge(full.ref, tlmps) %>% run.umap.10.20()
DimPlot(big.ref, reduction = "umap20", label = T, group.by = "sample.name") + NoLegend()
DimPlot(big.ref, reduction = "umap20", label = T, group.by = "combined.celltype") + NoLegend()
DimPlot(big.ref, reduction = "umap20", label = T, group.by = "TLMP_RNA_snn_res.1")

big.ref@meta.data$cell.name = colnames(big.ref)
big.ref@meta.data = big.ref@meta.data %>% left_join(all.dns@meta.data %>% select(cell.name, "all.dns.clusters"= RNA_snn_res.0.3), by = "cell.name")
rownames(big.ref@meta.data) = colnames(big.ref)

DimPlot(big.ref, reduction = "umap20", label = T, group.by = "all.dns.clusters")

#save figs

DimPlot(all.dns, reduction = "umap20", label = T)

FeaturePlot(all.dns, c("CD1A", "CD34", "CD33"), reduction = "umap20")

early.t.panel = c("CD1A", "CD2", "CD3D", "CD3E", "CD4", "CD7", "CD33", "CD34", "PTPRC", "KIT", "FLT3", "IL2RB", "IL7R", "CD5", "NFE2", "IGLL1")

early.t.panel[!early.t.panel %in% rownames(all.dns)]

directory.to.save = "/mnt/isilon/tan_lab/xuj5/HealthyReference/RScripts/07_TLMP_annotate_and_combine/"
p1=FeaturePlot(all.dns, early.t.panel, reduction = "umap20")

ggsave(plot = p1, filename = "early.dn.earlyT.panel.png", 
       path = directory.to.save, 
       device = "png", height = 15, width = 15)

ggsave(DimPlot(all.dns, reduction = "umap20", label = T),
       filename = "early.dn.clusters.png", 
       path = directory.to.save, 
       device = "png", height = 15, width = 15)

ggsave(VlnPlot(all.dns, early.t.panel, pt.size = 0.1),
       filename = "early.dn.vln.png", 
       path = directory.to.save, 
       device = "png", height = 15, width = 15)


ggsave(FeatureScatter(all.dns, "CD33", "CD34"),
       filename = "cd34_cd33.png", 
       path = directory.to.save, 
       device = "png", height = 15, width = 15)

ggsave(DimPlot(all.dns, reduction = "umap20", label = T, group.by = "seurat_clusters"),
       filename = "all.dn.clusters.png", 
       path = directory.to.save, 
       device = "png", height = 15, width = 15)

pdf(paste0(directory.to.save, "TLMP_sort_combined_withDN_andRef.pdf"))
DimPlot(tlmps ,label = T)
DimPlot(all.dns, reduction = "umap20", label = F, group.by = "sample.name")
DimPlot(all.dns, reduction = "umap20", label = T, group.by = "combined.celltype")
DimPlot(all.dns, reduction = "umap20", label = T, group.by = "TLMP_RNA_snn_res.1")
DimPlot(all.dns, reduction = "umap20", label = T, group.by = "Phase")
FeaturePlot(all.dns, reduction = "umap20", "G2M.Score")
DimPlot(big.ref, reduction = "umap20", label = T, group.by = "sample.name") + NoLegend()
DimPlot(big.ref, reduction = "umap20", label = T, group.by = "combined.celltype") + NoLegend()
DimPlot(big.ref, reduction = "umap20", label = T, group.by = "TLMP_RNA_snn_res.1")
DimPlot(big.ref, reduction = "umap20", label = T, group.by = "all.dns.clusters")
dev.off()

plot.degs = function(object, folder, object.name){
  
  directory.to.save = paste0("/mnt/isilon/tan_lab/xuj5/HealthyReference/RScripts/07_TLMP_annotate_and_combine/", folder, "/")
  dir.create(directory.to.save)
  
  tlmp.cluster.markers = FindAllMarkers(object, slot = "counts",
                                        logfc.threshold = 0.5, 
                                        latent.vars = "nCount_RNA", 
                                        return.thresh = 0.05)
  plot.list = list()
  for (i in 1:length(unique(tlmp.cluster.markers$cluster))){
    
    cluster.to.plot = unique((tlmp.cluster.markers$cluster))[i] %>% print()
    
    #filter
    cluster.makers = tlmp.cluster.markers[tlmp.cluster.markers$cluster == cluster.to.plot,]
    
    #volacnoPlot
    plot.list[[i]] = EnhancedVolcano(cluster.makers,
                                     lab = cluster.makers$gene,
                                     x = 'avg_logFC',
                                     y = 'p_val_adj', pCutoff = 0.05,
                                     FCcutoff = 0.5) + ggtitle(paste(object.name, "cluster", cluster.to.plot))
    
  }
  
  plot.list[[i+1]] = DimPlot(object, label = T, group.by = "seurat_clusters")
  
  ggsave(plot = cowplot::plot_grid(plotlist = plot.list),filename = paste0(object.name, "_cluster_DEGs.png"),device = "png", width = 20, height = 20, path = directory.to.save)
  
  ggsave(DimPlot(object, reduction = "umap20", label = T, group.by = "seurat_clusters"),
         filename = paste0(object.name, "_clustersonUMAP.png"), 
         path = directory.to.save, 
         device = "png", height = 15, width = 15)
}


plot.degs(tlmps, folder = "tlmp_cluster_markers", object.name = "tlmps only")
plot.degs(all.dns, folder = "all_dn_cluster_markers", object.name = "all dns")

#reassign clusters to big.ref
big.ref$updated.ids_5_14_21 = big.ref$combined.celltype
big.ref$updated.ids_5_14_21[big.ref$all.dns.clusters %in% c(0,1)] = "Pro-T"
big.ref$updated.ids_5_14_21[big.ref$all.dns.clusters %in% c(3,4)] = "Pre-T (Cycling)"
big.ref$updated.ids_5_14_21[big.ref$all.dns.clusters %in% c(2,5)] = "Pre-T"
big.ref$updated.ids_5_14_21[big.ref$all.dns.clusters %in% c(6)] = "TLMP"
big.ref$updated.ids_5_14_21[big.ref$all.dns.clusters %in% c(7)] = "pDC"

ggsave(DimPlot(big.ref, reduction = "umap20", label = T, group.by = "updated.ids_5_14_21") + ggtitle("full reference, 2000 VEGs, default settings"),
       filename = paste0("bigref_updated_ids_5_14_21.png"), 
       path = directory.to.save, 
       device = "png", height = 15, width = 15) 


saveRDS(big.ref, "/mnt/isilon/tan_lab/xuj5/HealthyReference/other_r_objects/thymus_bm/2021_5_14_full_ref_2000_defaultVEG")


