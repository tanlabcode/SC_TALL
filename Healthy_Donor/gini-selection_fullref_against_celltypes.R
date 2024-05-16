#gini selection
source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/RScripts/scDataAnalysis_Utilities.R')


#full.ref= readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/other_r_objects/tan_lab_combined_ref/4_23_full_ref_DN_early_late.RDS")
full.ref = readRDS("/mnt/isilon/tan_lab/xuj5/HealthyReference/other_r_objects/thymus_bm/2021_5_14_full_ref_2000_defaultVEG")

#must update column of metadata of updated celltypes
do.gini = function(seurat.rna, 
                   home.directory.to.save = paste0("/mnt/isilon/tan_lab/xuj5/HealthyReference/RScripts/08_full_ref_doGini/", Sys.Date(), "_GiniSelection_4000genes_start/")){
  
  source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/RScripts/scDataAnalysis_Utilities.R')
  
  #store orgiinal genes
  
  seurat.rna = seurat.rna %>% FindVariableFeatures(nfeatures = 4000) %>% ScaleData() %>% RunPCA(npcs = 30) %>% RunUMAP(dims = 1:20)
  
  original.vegs = VariableFeatures(seurat.rna)
  
  #save Gini Genes directory
  dir.create(home.directory.to.save)
  gini.genes.dir = paste0(home.directory.to.save, "2ndRoundGiniGenes&GO/") %>% print()
  dir.create(gini.genes.dir)

  
  print('getting og genes')
  #original genes
  original.vegs %>% as.data.frame() %>% write_delim(path = paste0(gini.genes.dir, "originalVEGs_", length(VariableFeatures(seurat.rna)), ".tsv"))
  gp1 = LabelPoints(plot=VariableFeaturePlot(seurat.rna), points=head(VariableFeatures(seurat.rna), 40), repel = T)  + 
    ggtitle(paste0("Original VEGs:", length(VariableFeatures(seurat.rna)), "genes"))
  ggsave(gp1, filename = paste0("originalVEGs_", length(VariableFeatures(seurat.rna)), ".png"),
         path = gini.genes.dir, height = 10, width =10)
  
  
  print('finding communities 1')
  #find communities
  #seurat.rna = seurat.rna %>% RunPCA(npcs = 50) %>% FindNeighbors(dims = 1:50, reduction = "pca") %>% FindClusters(res = 0.1) #find clusters
  #DimPlot(seurat.rna, label = T) + NoLegend()
  
  print('saving input 1')
  ggsave(filename = "input_to_first_gini.png", path = home.directory.to.save, plot = DimPlot(seurat.rna, label = T) + ggtitle("input to first Gini"))
  
  ## First gini selection
  print('doing 1st gini selection')
  gini_genes1 <- ifg_select(seurat.rna@assays$RNA@counts, 
                            seurat.rna$updated.ids_5_14_21, 
                            gini_cut_qt = 0.75, 
                            filePath = home.directory.to.save,
                            do_go = F)$include_g
  
  #set VEGs
  print('setting vegs, scaling data')
  VariableFeatures(seurat.rna) <- gini_genes1 
  
  seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                          vars.to.regress = c('perc.mito', 'nCount_RNA'))
  seurat.rna <- RunPCA(seurat.rna, npc = 30, verbose = F, 
                       features = VariableFeatures(seurat.rna))
  print('runing umap')
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:20, reduction.name = "umap_gini_1")
  #seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)

  
  p1 =  DimPlot(seurat.rna, label = T, reduction = "umap_gini_1")+ ggtitle("communities after gini 1st round")
  p2 =  DimPlot(seurat.rna, group.by = "updated.ids_5_14_21", label = T, reduction = "umap_gini_1") + ggtitle("known cell types after gini 1st round") + NoLegend()
  
  print('saving first gini output')
  ggsave(cowplot::plot_grid(plotlist = list(p1, p2)),
         filename = paste0("AfterGini1stRound_UMAP20.png"), 
         path = home.directory.to.save, 
         device = "png", height = 10, width = 15)
  
  #1st gini
  
  print('storing gini genes')
  gini_genes1 %>% as.data.frame() %>% write_delim(path = paste0(gini.genes.dir, "1stGiniVEGs_", length(gini_genes1), ".tsv"))
  gp2 = LabelPoints(plot=VariableFeaturePlot(seurat.rna), points=head(gini_genes1, 40), repel = T) + ggtitle(paste0("1st Gini Selection:", length(gini_genes1), "genes"))
  ggsave(gp2, filename = paste0("1stGiniVEGs_", length(gini_genes1), ".png"),
         path = gini.genes.dir, height = 10, width =10)
  
  
  
  
  ## Second gini selection
  #find communities
  print('building graph for gini 2')
  #seurat.rna <- FindNeighbors(seurat.rna, dims = 1:20, reduction = 'pca')
  #seurat.rna <- FindClusters(seurat.rna, res = 0.1)
  
  #find gini genes
  gini_genes2 <- ifg_select(seurat.rna@assays$RNA@counts, 
                           seurat.rna$updated.ids_5_14_21, 
                           gini_cut_qt = 0.75,
                           filePath = home.directory.to.save)$include_g
  #set gini genes
  print('setting vegs for gini 2')
  VariableFeatures(seurat.rna) <- gini_genes2
  
  #scale Data
  print('scaling data')
  seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                          vars.to.regress = c('perc.mito', 'nCount_RNA'))
  
  #Run PCA
  print('ruunning pca and umap')
  seurat.rna <- RunPCA(seurat.rna, npc = 50, verbose = F, 
                       features = VariableFeatures(seurat.rna))
  
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:50, reduction.name = "umap_gini_2")
  #seurat.rna <- RunTSNE(seurat.rna, dims = 1:50)
  seurat.rna <- FindNeighbors(seurat.rna, dims = 1:50, reduction = 'pca')
  seurat.rna <- FindClusters(seurat.rna, res = 0.6)
  
  
  p3 =   DimPlot(seurat.rna, group.by = "updated.ids_5_14_21", reduction = "umap_gini_2", label = T) + NoLegend() + ggtitle("celltypes after 2nd Gini (UMAP50)")
  p4 =   DimPlot(seurat.rna, reduction = "umap_gini_2", label = T) + ggtitle("communities after second gini (UMAP50)")
  print('second gini output saving graph')
  ggsave(cowplot::plot_grid(plotlist = list(p3, p4)),
         filename = paste0("AfterGini2ndRound_UMAP50.png"), 
         path = home.directory.to.save, 
         device = "png", height = 10, width = 15)
  
  
  directory.to.save = paste0(home.directory.to.save, "2ndGiniExploreUMAPinput/")
  dir.create(directory.to.save)
  
  print('making umaps over many PCs')
  npc.input = c(10, 15, 20, 30, 40, 50)
  for (i in 1:length(npc.input)){
    
    reduction.name = paste0("umap_gini_2_", npc.input[i]) %>% print()
    
    seurat.rna <- RunUMAP(seurat.rna, dims = 1:npc.input[i], reduction.name = reduction.name)
    
    
    png(paste0(directory.to.save, "combinedref-secondGini-by-population_UMAP_", npc.input[i], ".png"))
    print(DimPlot(seurat.rna, group.by = "updated.ids_5_14_21", reduction = reduction.name, label = T) + NoLegend())
    dev.off()
    
    
    png(paste0(directory.to.save, "combinedref-secondGini_UMAP", npc.input[i], ".png"))
    print(DimPlot(seurat.rna, reduction = reduction.name, label = T))
    dev.off()
    
  }
  
  #make plot of all
  print('saving all umaps in plotgrid')
  plot.list = list()
  for (i in 1:length(npc.input)){
    
    reduction.name = paste0("umap_gini_2_", npc.input[i]) %>% print()
    
    plot.list[[i]] = (DimPlot(seurat.rna, group.by = "updated.ids_5_14_21", reduction = reduction.name, label = T) + NoLegend()) + ggtitle(reduction.name)
    
  }
  plot.list[[i+1]]= ElbowPlot(seurat.rna, ndims = 50)
  
  ggsave(cowplot::plot_grid(plotlist = plot.list), filename = "after2nd_umap10-40.png", height = 15, width = 15, path = directory.to.save)
  

  #2nd gini veg plot
  print('second gini gini genes saving tsv')
  gini_genes2 %>% as.data.frame() %>% write_delim(path = paste0(gini.genes.dir, "2ndGiniVEGs_", length(gini_genes2), ".tsv"))
  gp3 = LabelPoints(plot=VariableFeaturePlot(seurat.rna), points=head(gini_genes2, 40), repel = T) + ggtitle(paste0("2nd Gini Selection:", length(gini_genes2), "genes"))
  ggsave(gp3, filename = paste0("2ndGiniVEGs_", length(gini_genes2), ".png"),
         path = gini.genes.dir, height = 10, width =10)
  
  #allrounds [plot grid]
  print('saving veg plot for all gini genes')
  ggsave(cowplot::plot_grid(plotlist = list(gp1, gp2, gp3)), filename = paste0("allRoundsGenes.png"),
         path = gini.genes.dir, height = 10, width =20)
  
  print('doing go')
  go.anal.og = do_GO(original.vegs, type = "BP", qCutoff = 0.05,
                     organism = "hsa")
  go.anal.gini.2 = do_GO(gini_genes2, type = "BP", qCutoff = 0.05,
                         organism = "hsa")
  
  
  
  thrown.out = original.vegs[!original.vegs %in% gini_genes2] 
  
  thrown.out.go = do_GO(thrown.out, type = "BP", qCutoff = 0.05,
                        organism = "hsa")
  
  thrown.out %>% as.data.frame() %>% write_delim(path = paste0(gini.genes.dir, "thrownout2_", length(thrown.out), ".tsv"))
  
  print('writing GO')
  write_tsv(go.anal.og@result,  paste0(gini.genes.dir, "original_enriched_GO.tsv"))
  write_tsv(go.anal.gini.2@result, paste0(gini.genes.dir, "2ndGini_enriched_GO.tsv"))
  write_tsv(thrown.out.go@result, paste0(gini.genes.dir, "thrownout2_enriched_GO.tsv"))
  
  return(seurat.rna)
}

ref.after.gini = do.gini(seurat.rna = full.ref)

VariableFeatures(ref.after.gini) %>% length()


saveRDS(ref.after.gini, "/mnt/isilon/tan_lab/xuj5/HealthyReference/other_r_objects/thymus_bm/2021_5_14_full_ref_1048_GiniVEG")
