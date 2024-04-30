#import final matrixes

setwd("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/01_mkSeurats")

source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/RScripts/scDataAnalysis_Utilities.R')
library(Seurat)
library(RColorBrewer)
library(tidyverse)
library(data.table)

sample.df <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/sampledfs/merged_dfs/Apr9_65samples_CITESeq_sampledf_with_metadata_merged.rds")

sample.df

files = list.files('/mnt/isilon/tan_lab/xuj5/ETP_ALL/SeuratObjects', full.names = T) %>% print()
length(files)

make.big.seurat = function(paths){
  print(paste("reading in", paths[1], "which is", 1, "out of", length(paths)))
  seurat.big = readRDS(paths[1]) #initilize first seurat
  
  for(i in 2:length(paths)){
    
    print(paste("reading in", paths[i], "which is", i, "out of", length(paths)))
    
    to.add =readRDS(paths[i]) #read in next seurat
    
    print("merging objects")
    seurat.big = merge(seurat.big, to.add)
      
  }
  
  return(seurat.big)
  
}

all.etp.65 = make.big.seurat(files)

make.global.umap = function(seurat.obj){
  
  print("finding var features")
  seurat.obj <- FindVariableFeatures(seurat.obj, nfeatures = 5000)
  
  #find top 5000 features
  print("taking features with expression in >0.01% cells ")
  variable.features <- seurat.obj@assays$RNA@counts[VariableFeatures(seurat.obj),]
  
  perc.cells.express = rowSums(variable.features >= 1) /ncol(variable.features)
  
  table(perc.cells.express > 0.01)
  
  VariableFeatures(seurat.obj) = names(perc.cells.express)[perc.cells.express > 0.01] #take only those expressed in at least 1% of cells
  
  print(paste("retained", VariableFeatures(seurat.obj) %>% length(), "features"))
  
  print('running PCA + UMAP')
  seurat.obj <- seurat.obj %>% 
    ScaleData(features = VariableFeatures(seurat.obj), vars.to.regress = c("perc.mito","nCount_RNA")) %>% 
    RunPCA(npcs = 100) %>% RunUMAP(dims = 1:50, reduction.name = "umap_noCC_HS_Reg")
  
  
  
  #get heatshock and cell cycle score
  print("cellcycle regression")
  #cell cycle regression
  cycle3 = fread('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/Signature_GeneList/regev_lab_cell_cycle_genes.txt', header = F)$V1
  s.genes = cycle3[1:43]
  g2m.genes = cycle3[44:97]
  seurat.obj = CellCycleScoring(object = seurat.obj, s.features = s.genes,
                                g2m.features = g2m.genes, set.ident = TRUE)
  
  #heat shock regression
  heat_shock_gene = fread('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/Signature_GeneList/heat_shock_geneList.txt')
  heat_shock_gene = heat_shock_gene$`Approved symbol`
  seurat.obj = AddModuleScore(seurat.obj, features = list(heat_shock_gene),
                              name = 'HeatShock.Score')
  
  print("rerunning UMAP")
  seurat.obj <- seurat.obj %>% 
    ScaleData(features = VariableFeatures(seurat.obj), vars.to.regress = c("S.Score", "G2M.Score", "perc.mito","HeatShock.Score1", "nCount_RNA")) %>% 
    RunPCA(npcs = 100)%>% RunUMAP(dims = 1:50, reduction.name = "umap_with_CC_HS_Reg")
  
  print("returning seurat obj")

  
  #save some figures
  
  print('saving figures')
  directory.to.save = paste0("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/02_combine_objects_mkUMAP/", length(seurat.obj$sample.name %>% unique()),
                             "_samples_UMAP")
  
  dir.create(directory.to.save)
  
  p1 = DimPlot(seurat.obj, group.by = "sample.name", label = T, reduction = "umap_noCC_HS_Reg") + NoLegend() + ggtitle("no cell cycle / heat shock reg")
  p2 = DimPlot(seurat.obj, group.by = "sample.name", label = T, reduction = "umap_with_CC_HS_Reg") + NoLegend() + ggtitle("with cell cycle / heat shock reg")
  p3 = DimPlot(seurat.obj, group.by = "sample.group", label = T, reduction = "umap_noCC_HS_Reg")  + ggtitle("no cell cycle / heat shock reg")
  p3 = DimPlot(seurat.obj, group.by = "sample.group", label = T, reduction = "umap_with_CC_HS_Reg") + ggtitle("with cell cycle / heat shock reg")
  
  ggsave(plot = cowplot::plot_grid(plotlist = list(p1, p2, p3, p4)),filename = "UMAP50_grid.png",device = "png", width = 20, height = 20, path = directory.to.save)

  return(seurat.obj)
}




all.etp.65= make.global.umap(all.etp.65)

save.figures = function(seurat.obj){
  
  
  print('saving figures')
  directory.to.save = paste0("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/02_combine_objects_mkUMAP/", length(seurat.obj$sample.name %>% unique()),
                             "_samples_UMAP")
  
  dir.create(directory.to.save)
  
  p1 = DimPlot(seurat.obj, group.by = "sample.name", label = T, reduction = "umap_noCC_HS_Reg") + NoLegend() + ggtitle("no cell cycle / heat shock reg")
  p2 = DimPlot(seurat.obj, group.by = "sample.name", label = T, reduction = "umap_with_CC_HS_Reg") + NoLegend() + ggtitle("with cell cycle / heat shock reg")
  p3 = DimPlot(seurat.obj, group.by = "sample.group", label = T, reduction = "umap_noCC_HS_Reg")  + ggtitle("no cell cycle / heat shock reg")
  p4 = DimPlot(seurat.obj, group.by = "sample.group", label = T, reduction = "umap_with_CC_HS_Reg") + ggtitle("with cell cycle / heat shock reg")
  
  print("ggsave")
  ggsave(plot = cowplot::plot_grid(plotlist = list(p1, p2, p3, p4)),filename = "UMAP50_grid.png",device = "png", width = 20, height = 20, path = directory.to.save)
  
}

save.figures(all.etp.65)

all.etp.65$sample.group[all.etp.65$sample.group == "ETP High MRD no relapse "] = "ETP High MRD no relapse"
all.etp.65$sample.group %>% table()

sample.group.plot = tibble(sample.group = unique(all.etp.65$sample.group),
                     sample.group.short = c("AML", "HealthyThymus", "MPAL", "ETP", "ETP", "ETP", "Non-ETP"))

all.etp.65@meta.data$cell.name = rownames(all.etp.65@meta.data)
all.etp.65@meta.data = all.etp.65@meta.data %>% left_join(sample.group.plot)
rownames(all.etp.65@meta.data)=all.etp.65$cell.name

directory.to.save = paste0("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/02_combine_objects_mkUMAP/", length(all.etp.65$sample.name %>% unique()),
                           "_samples_UMAP")

graph.composition = all.etp.65@meta.data %>% group_by(sample.group.short) %>% summarize(n = length(unique(sample.name))) %>% mutate(z = paste(n, sample.group.short))

title = paste(ncol(all.etp.65), "cells from", unique(all.etp.65$sample.name) %>% length(), "samples:",
              graph.composition$z %>% paste(collapse = "; "))

p1 = DimPlot(all.etp.65, group.by = "sample.name", label = T, reduction = "umap_noCC_HS_Reg") + NoLegend() + labs(title = title, subtitle = "no cell cycle / heat shock reg")+ NoAxes()
p2 = DimPlot(all.etp.65, group.by = "sample.name", label = T, reduction = "umap_with_CC_HS_Reg") + NoLegend() + labs(title = title, subtitle ="with cell cycle / heat shock reg")+ NoAxes()
p3 = DimPlot(all.etp.65, group.by = "sample.group.short", label = F, reduction = "umap_noCC_HS_Reg")  + labs(title = title, subtitle ="no cell cycle / heat shock reg") +
  theme(legend.position = "bottom") + NoAxes()
p4 = DimPlot(all.etp.65, group.by = "sample.group.short", label = F, reduction = "umap_with_CC_HS_Reg") + labs(title = title, subtitle = "with cell cycle / heat shock reg") +
  theme(legend.position = "bottom") + NoAxes()

ggsave(plot = cowplot::plot_grid(plotlist = list(p1, p2, p3, p4)),filename = "UMAP50_grid.png",device = "png", width = 20, height = 20, path = directory.to.save)

saveRDS(all.etp.65, "/mnt/isilon/tan_lab/xuj5/ETP_ALL/other_r_objects/combined_seurat_objects/May5_2021_65_samples.rds")
