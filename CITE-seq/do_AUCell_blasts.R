library(AUCell)
dir.save = paste0(folders$final.r.scripts, "/29_doAUCell/")
dir.create(dir.save)
all.etp = readRDS(files$all.etp.with.gini.annos)


marker.list <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/28_AUCellGeneSets/cell.type.short_11603_cells/all_markers_list.rds")

all.markers.all.celltypes <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/28_AUCellGeneSets/cell.type.short.against.2.neighbors/all.markers.all.celltypes.rds")

get.downsample.matrix.all.etp = function(seurat.obj, cell.id.column, amount.keep = NA){
  
  downsample.matrix = seurat.obj@meta.data %>% group_by(.data[[cell.id.column]]) %>% summarize(n = length(nCount_RNA))
  downsample.matrix = downsample.matrix %>% left_join(seurat.obj@meta.data %>% dplyr::select(.data[[cell.id.column]]) %>% unique())
  downsample.matrix$amount.to.keep = amount.keep
  return(downsample.matrix)
  }



downsample.matrix = get.downsample.matrix.all.etp(all.etp, cell.id.column = "sample.name", amount.keep = 1000)

get.keepers.all.etp = function(downsample.matrix, seurat.obj){
  
  to.keep.global = tibble(cell.name = NULL,
                          sample.name = NULL)
  
  for (i in 1:length(unique(downsample.matrix$sample.name))){
    
    print(i)
    sample = unique(downsample.matrix$sample.name)[i] %>% print()
    number.in.ref = downsample.matrix$n[i] %>% print()
    number.to.keep = downsample.matrix$amount.to.keep[i] %>% print()
    
    
    cells.of.cell.type = seurat.obj$cell.name[seurat.obj$sample.name == sample]
    
    if(number.in.ref >number.to.keep){
      
      keepers = sample(cells.of.cell.type, size = number.to.keep, replace = F)
      
    }else{
      
      keepers = cells.of.cell.type
      
    }
    
    keepers.table  = tibble(cell.name = keepers,
                            sample.name = sample)
    
    to.keep.global = rbind(to.keep.global, keepers.table)
    
  }
  
  to.keep.global$sample.name %>% table() %>% print()
  
  return(to.keep.global)
  
}


keepers = get.keepers.all.etp(downsample.matrix, all.etp)

downsample.etp = subset(all.etp, cell.name %in% keepers$cell.name)

cells_rankings_etp <- AUCell_buildRankings(downsample.etp@assays$RNA@counts, nCores = 12)

au.cell.list = find.all.markers.to.aucell.list(find.all.markers.list = marker.list$house.genes.removed)

cells_AUC <- AUCell_calcAUC(au.cell.list, cells_rankings_etp, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC_TFs <- AUCell_calcAUC(find.all.markers.to.aucell.list(find.all.markers.list = marker.list$tfs), cells_rankings_etp, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC_nn <- AUCell_calcAUC(nn.markers.list.to.au.cell.list(all.markers.all.celltypes), cells_rankings_etp, aucMaxRank=nrow(cells_rankings)*0.05)

downsample.etp[["AUCell"]] <- CreateAssayObject(cells_AUC@assays@data$AUC)
downsample.etp[["AUCell_TF"]] <- CreateAssayObject(cells_AUC_TFs@assays@data$AUC)
downsample.etp[["AUCell_nn"]] <- CreateAssayObject(cells_AUC_nn@assays@data$AUC)

DefaultAssay(downsample.etp)= "AUCell"

downsample.etp = downsample.etp %>% ScaleData()


DoHeatmap(toplot, features = rownames(toplot@assays$AUCell), slot = "scale.data", label = T, raster = T, group.by = "level.1.anno.viscello",
          cells = sample(colnames(toplot), size = 0.10 * ncol(toplot), replace = F)) + 
  scale_fill_gradient(low = heat.colors(n=10)[10], high = heat.colors(n=10)[2])


all.etp.blasts = subset(toplot, level.1.anno.viscello == "Blast")

DoHeatmap(toplot.blasts, features = rownames(toplot.blasts@assays$AUCell), slot = "scale.data", label = T, raster = T, group.by = "predicted.id",
          cells = sample(colnames(toplot.blasts), size = 0.10 * ncol(toplot.blasts), replace = F)) + 
  scale_fill_gradient(low = heat.colors(n=10)[10], high = heat.colors(n=10)[2])

