#projection function

#noccreg
library(Seurat)
library(tidyverse)
source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/RScripts/scDataAnalysis_Utilities.R')
source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/JX_scFunctions.R')


project.samples.v4 = function(big.seurat = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/other_r_objects/combined_seurat_objects/2021_7_5_60_leukemia_samples_with_level1anno.rds"), 
                           reference.traj = readRDS("/mnt/isilon/tan_lab/xuj5/HealthyReference/other_r_objects/thymus_bm/2021_6_14_full.ref.full.genes.update.etp.rds"),
                           vegs = NA,
                           npc = 30,
                           min.dist = 0.5,
                           directory.to.save =  paste0("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/14_project_using_seurat_v4/", Sys.Date(), "_ProjectionFiles_", length(big.seurat$sample.name %>% unique()), "_Samples_", npc, "_PCs_", length(VariableFeatures(reference.traj)), "_VEGs/")
){
  
  
  dir.create(directory.to.save)
  
  reference.traj <- reference.traj
  reference.traj = reference.traj %>% NormalizeData() %>% ScaleData(features = VariableFeatures(reference.traj), vars.to.regress = c("nCount_RNA", "perc.mito")) %>% 
    RunPCA(npcs = 100)%>% 
    RunUMAP(dims = 1:npc, return.model = T,min.dist = 0.5, reduction.name = "umap_4_projection", umap.method = "uwot", reduction = "pca", assay = "RNA")
  
  
  DimPlot(reference.traj, group.by = "cell.type.short", reduction = "umap_4_projection", label = T) + NoLegend() 
  ggsave(filename = "ref.traj.png", path = directory.to.save, device = "png", height = 10, width = 10)
  
  
  reference.traj <- FindNeighbors(
    object = reference.traj,
    reduction = "pca",
    dims = 1:npc,
    graph.name = "pca.annoy.neighbors", 
    k.param = 50,
    cache.index = TRUE,
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  saveRDS(reference.traj, paste0(directory.to.save, "ref.traj.umap4projection.rds"))
  
  #6 objects
  #setwd("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/06_check_level_one_anno/2021-06-01_level_1_anno_6_Samples_by_group/objects_by_group")
  
  all.etp = big.seurat
  
  #query.objects = list.files(full.names = T, path ="/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/06_check_level_one_anno/2021-06-01_level_1_anno_6_Samples_by_group/objects_by_group")
  
  query.objects = SplitObject(all.etp, split.by = "sample.group")
  
  
  all.sample.coords = NULL
  for (i in 1:length(query.objects)){
    
    group.to.proj= query.objects[[i]]
    
    print('splitting object')
    samples.to.project <- SplitObject(group.to.proj, split.by = "sample.name")
    
    print('normalizing data')
    samples.to.project <- lapply(X = samples.to.project, FUN = NormalizeData, verbose = FALSE)
    
    print("finding anchors")
    anchors <- list()
    for (i in 1:length(samples.to.project)) {
      print(paste("anchoring", i))
      anchors[[i]] <- FindTransferAnchors(
        reference = reference.traj,
        query = samples.to.project[[i]],
        k.filter = NA,
        reference.reduction = "pca", 
        reference.neighbors = "pca.annoy.neighbors", 
        dims = 1:npc
      )
    }
    
    
    for (i in 1:length(samples.to.project)) {
      print(paste("mapping", i))
      samples.to.project[[i]] <- MapQuery(
        anchorset = anchors[[i]], 
        query = samples.to.project[[i]],
        reference = reference.traj, 
        refdata = list(t.traj.ptime = "curve2",
                       m.traj.ptime = "curve5",
                       cell.type.short = "cell.type.short",
                       cell.type.binned = "cell.type.binned",
                       trajectory = "trajectory"),
        reference.reduction = "pca",
        reduction.model = "umap_4_projection"
      )
    }
    
    
    group.to.proj <- merge(samples.to.project[[1]], samples.to.project[2:length(samples.to.project)], merge.dr = c("ref.umap", "umap_noCC_HS_Reg", "umap_with_CC_HS_Reg"))
    
    sample.group= group.to.proj$sample.group %>% unique()
    
    #save object
    object.name = paste0(directory.to.save, sample.group[1], "_", length(group.to.proj$sample.name %>% unique()), "_Samples_", npc, "_PCs_", length(VariableFeatures(reference.traj)), "_VEGs_projection.rds")
    saveRDS(group.to.proj, file = object.name)
    
    sample.coords = group.to.proj@reductions$ref.umap@cell.embeddings %>% as_tibble(rownames = "cell.name")
    sample.coords = sample.coords %>% left_join(group.to.proj@meta.data, by = "cell.name")
    
    all.sample.coords = rbind(all.sample.coords, sample.coords)
    
    
  }
  
  
  #save all umap coords & predicted ids
  
  ref.coords = reference.traj@reductions$umap_4_projection@cell.embeddings %>% as_tibble(rownames = "cell.name")
  write_csv(ref.coords, file = paste0(directory.to.save, "ref.coords.csv"))
  
  
  write_csv(all.sample.coords, file = paste0(directory.to.save, "umap_coords_2_with_metadata.csv"))
  
  colnames(ref.coords) <- c("cell.name", "refUMAP_1", "refUMAP_2")
  ref.coords <- ref.coords %>%left_join(reference.traj@meta.data, by = "cell.name")
  
  write_csv(ref.coords, file = paste0(directory.to.save, "ref.coords.csv"))
  
  base.umap <- ggplot(ref.coords, aes(x=umap_4_projection_1, y = umap_4_projection_2)) + geom_point(color = "gray", alpha = 1, size = 0.5) + theme_minimal()
  
  umap.coords.table <- all.sample.coords %>% arrange(sample.group)
  
  umap.coords.table$sample.name <- factor(umap.coords.table$sample.name, levels = umap.coords.table$sample.name %>% unique())
  
  base.umap + geom_point(data = umap.coords.table, alpha = 1, aes(x= refUMAP_1, y = refUMAP_2, color = sample.group), size = 0.1) + facet_wrap(~sample.name) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(filename = "all.sample.proj.dot.png", path = directory.to.save, device = "png", height = 10, width = 10)
  
  base.umap + geom_density_2d(data = umap.coords.table, aes(x= refUMAP_1, y = refUMAP_2, color = sample.group), binwidth = 0.01) + facet_wrap(~sample.group) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(filename = "all.group.dens.png", path = directory.to.save, device = "png", height = 10, width = 10)
  
  base.umap + geom_density_2d(data = umap.coords.table, alpha = 1, aes(x= refUMAP_1, y = refUMAP_2, color = sample.group), binwidth = 0.01) + facet_wrap(~sample.name) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(filename = "all.sample.proj.dens.png", path = directory.to.save, device = "png", height = 10, width = 10)

}

project.samples.v4(reference.traj = readRDS(paste0(folders$final.r.scripts, "/22_doSlingshot", "/ref-traj-all-curve-pseudotime-by-umap.rds")))

make.plots = function(directory.to.save = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/14_project_using_seurat_v4/2021-06-24_ProjectionFiles_60_Samples_30_PCs_1037_VEGs/"){
  
  ref.coords = read_csv(paste0(directory.to.save, "ref.coords.csv"))
  all.sample.coords = read_csv(paste0(directory.to.save, "umap_coords_2_with_metadata.csv"))
  
  base.umap <- ggplot(ref.coords, aes(x=UMAP_1, y = UMAP_2)) + geom_point(color = "gray", alpha = 1, size = 0.5) + theme_minimal()
  
  umap.coords.table <- all.sample.coords %>% arrange(sample.group)
  
  umap.coords.table$sample.name <- factor(umap.coords.table$sample.name, levels = umap.coords.table$sample.name %>% unique())
  
  base.umap + geom_point(data = umap.coords.table, alpha = 1, aes(x= refUMAP_1, y = refUMAP_2, color = sample.group), size = 0.1) + facet_wrap(~sample.name) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(filename = "all.sample.proj.dot.png", path = directory.to.save, device = "png", height = 10, width = 10)
  
  base.umap + geom_density_2d(data = umap.coords.table, aes(x= refUMAP_1, y = refUMAP_2, color = sample.group), binwidth = 0.01) + facet_wrap(~sample.group) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(filename = "all.group.dens.png", path = directory.to.save, device = "png", height = 10, width = 10)
  
  base.umap + geom_density_2d(data = umap.coords.table, alpha = 1, aes(x= refUMAP_1, y = refUMAP_2, color = sample.group), binwidth = 0.01) + facet_wrap(~sample.name) +
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  ggsave(filename = "all.sample.proj.dens.png", path = directory.to.save, device = "png", height = 10, width = 10)
}



make.plots()



annotated.objects = list.files(directory.to.save, pattern ="projection.rds$") %>% print()

all.sample.coords = NULL
for (i in 1:length(annotated.objects)){
  print(i)
  group.to.proj = readRDS(annotated.objects[i])
  
  
  sample.coords = group.to.proj@reductions$ref.umap@cell.embeddings %>% as_tibble(rownames = "cell.name")
  sample.coords = sample.coords %>% left_join(group.to.proj@meta.data, by = "cell.name")
  
  all.sample.coords = rbind(all.sample.coords, sample.coords)
  
  
}


