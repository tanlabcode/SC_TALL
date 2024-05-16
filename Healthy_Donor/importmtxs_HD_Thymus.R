#process all HD into seurat objects
setwd("/mnt/isilon/tan_lab/xuj5/HealthyReference/RScripts")

source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/RScripts/scDataAnalysis_Utilities.R')
library(Seurat)
library(RColorBrewer)
library(tidyverse)
library(data.table)


print("done loading packages")

#set data.directories and dates----

sample.df = readRDS("/mnt/isilon/tan_lab/xuj5/HealthyReference/sample_dfs/45_samples_all_healthydonors.rds")
thymus.df = read_csv("/mnt/isilon/tan_lab/xuj5/HealthyReference/RScripts/02_process_raw_data_mkObject/thymus_samples.csv")


#create list of samples organized by date----

file.exists(thymus.df$matrix.location)


#read in data for sample i----

directories = list()
directories$seurat.object.directory <- '/mnt/isilon/tan_lab/xuj5/HealthyReference/SeuratObjects/Thymus/' %>% print()
directories$figure.directory <- paste0('/mnt/isilon/tan_lab/xuj5/HealthyReference/QCFigures') %>% print()
directories$figure.basic.qc <- paste0(directories$figure.directory, "/basic_qc/") %>% print()
directories$figure.doublet.plot <- paste0(directories$figure.directory, "/doublet_plot/") %>% print()
directories$figure.umap.clusters <- paste0(directories$figure.directory, "/umap50_clusters/") %>% print()
directories$figure.umap.gene.features <- paste0(directories$figure.directory, "/umap50_gene_features/") %>% print()
directories$figure.umap.adt.features <- paste0(directories$figure.directory, "/umap50_adt_features/") %>% print()
directories$figure.adt.violin.plot <-paste0(directories$figure.directory, "/adt_violin_plot/") %>% print()
directories$figure.all.features <- paste0(directories$figure.directory, "/umap50_adt_gene_features/") %>% print()
directories$figure.combined.pdf <- paste0(directories$figure.directory, "/combinedpdfs_QC_UMAP/") %>% print()
directories$filter.stats.directory <- paste0(directories$figure.directory, "/filterstats/") %>% print()
directories$all.raw.data.directory <- '/mnt/isilon/tan_lab/xuj5/HealthyReference/rawdata/'
directories$sample.df.directory <- '/mnt/isilon/tan_lab/xuj5/HealthyReference/sampledfs/' %>% print()

lapply(X = directories, dir.create)

sample.df = read_csv("/mnt/isilon/tan_lab/xuj5/HealthyReference/RScripts/02_process_raw_data_mkObject/thymus_samples.csv")

sample.df$starting.cells = sample.df$starting.genes = 
  sample.df$final.cells = sample.df$hb.filter = sample.df$mito.filter = 
  sample.df$nCount_1500.filter = sample.df$doublets = sample.df$nCount40000.filter = sample.df$median.nCount =sample.df$filter.adt.2000.hypothetical =
  NA

#make seurat obj for all leukemia objects
for (counter in 1:6){
  if(sample.df$CITEseq[counter]==T){
    
    data <-  Read10X(data.dir = sample.df$matrix.location[counter])
    sampleID <- sample.df$sample.name[counter] %>% print()
    sample.group <- sample.df$group[counter] %>% print()
    sequence.date <- sample.df$sequence.date[counter]
    process.date <- sample.df$process.date[counter]
    
    print(paste("processing sample", sampleID, "which is", counter, "out of", length(sample.df$sample.name), "samples"))
    
    
    colnames(data$`Gene Expression`) = paste0(sampleID, '_', colnames(data$`Gene Expression`))
    colnames(data$`Antibody Capture`) = paste0(sampleID, '_', colnames(data$`Antibody Capture`))
    
    original.dim <- dim(data$`Gene Expression`) %>% print()
    
    sample.df$starting.cells[counter] = ncol(data$`Gene Expression`)
    sample.df$starting.genes[counter] = nrow(data$`Gene Expression`)
    sample.df$hb.filter[counter] = sum(data$`Gene Expression`['HBB', ] > 3)
    
    print("starting to filter data")
    
    #filter hemoglobin----
    if('HBB' %in% rownames(data$`Gene Expression`)){
      hbb_exp = data$`Gene Expression`['HBB', ]
      data$`Gene Expression` = data$`Gene Expression`[, hbb_exp < 3]
    }
    
    if(sum(hbb_exp>3)>0){
      filter.rbc <- table(hbb_exp <3)
      print(paste("filtered", filter.rbc[1], "rbcs out, kept", filter.rbc[2], "cells"))
    }else{
      filter.rbc <- c(0, table(hbb_exp <3))
      print(paste("filtered", filter.rbc[1], "rbcs out, kept", filter.rbc[2], "cells"))
    }
    #filter by count between 1500 and 40,000
    nCount = Matrix::colSums(data$`Gene Expression`) #totalreads per cell
    
    sample.df$nCount_1500.filter[counter] = sum(nCount<1000)
    sample.df$nCount40000.filter[counter] = sum(nCount>40000)
    
    
    data$`Gene Expression` = data$`Gene Expression`[, nCount < 40000 & nCount > 1500]
    
    sample.df$median.nCount[counter] = median(Matrix::colSums(data$`Gene Expression`))
    
    filter.nCount <- table(nCount >1500 & nCount <40000) 
    print(paste("filtered", filter.nCount[1], "cells out, kept", filter.nCount[2], "cells"))
    
    #filter by perc.mito > 0.1
    mito.features <- grep(pattern = "^MT-", 
                          x = rownames(x = data$`Gene Expression`), value = TRUE)
    
    perc.mito = Matrix::colSums(data$`Gene Expression`[mito.features, ])/Matrix::colSums(data$`Gene Expression`)
    
    #plot all mito features
    #par(mfrow=c(3,5))
    #for (i in 1:length(mito.features)){
    #data$`Gene Expression`[mito.features[i], ] %>% hist(main = mito.features[i], xlim = c(1,500))
    #}
    
    #hist(perc.mito, main = sampleID) 
    
    #quantile(perc.mito)
    
    filter.mito <- table(perc.mito <= 0.1) %>% print()
    
    sample.df$mito.filter[counter] = sum(perc.mito >= 0.1)
    
    data$`Gene Expression` = data$`Gene Expression`[, perc.mito <= 0.1]
    perc.mito = perc.mito[perc.mito <= 0.1]
    
    #Subset Antibody Capture cells
    data$`Antibody Capture` = data$`Antibody Capture`[,colnames(data$`Antibody Capture`) %in% colnames(data$`Gene Expression`)]
    
    #filter by ncount_adt < 2000
    nCount_ADT <- Matrix::colSums(data$`Antibody Capture`)
    #data$`Antibody Capture` = data$`Antibody Capture`[, nCount_ADT < 2000]
    
    filter.adt <- table(nCount_ADT < 2000)
    sample.df$filter.adt.2000.hypothetical[counter] = sum(nCount_ADT<2000)
    #Subset Gene by Antibody counts
    
    data$`Gene Expression` = data$`Gene Expression`[,colnames(data$`Gene Expression`) %in% colnames(data$`Antibody Capture`)]
    
    #final filter stats
    filter.stats <- c(original.dim, filter.rbc, filter.nCount, filter.mito, filter.adt) %>% 
      matrix(ncol = 2, byrow = T, dimnames = list(c(paste(sampleID,"original.dim"), "filter.rbc", "filter.nCount", "filter.mito", "filter.adt(hypothetical)"), c("filtered", "kept"))) %>% print()
    
    print("writing filter statistics")
    write.csv(filter.stats, file = paste0(directories$filter.stats.directory, sampleID, "_filterstats"))
    
    
    
    #Create Seurat Object
    seurat.ETP = CreateSeuratObject(counts = data$`Gene Expression`) #first gene exp
    seurat.ETP[['ADT']] = CreateAssayObject(counts = data$`Antibody Capture`) #then ADT
    
    cnames = colnames(seurat.ETP) #cell names
    tmp.mito = data.table('perc' = perc.mito, 'cname' = names(perc.mito)) #datatable, perc mito and cells as column names
    setkey(tmp.mito, cname) #fixes key
    seurat.ETP@meta.data[['perc.mito']] = tmp.mito[cnames]$perc #set seurat object metadata for perc.mito
    
    #create metadata
    print("writing sample.name, sample.group, sequence.date, process.date, to sample metadata")
    seurat.ETP@meta.data[['sample.name']] = sampleID
    seurat.ETP@meta.data[['sample.group']] = sample.group
    seurat.ETP@meta.data[['sequence.date']] = sequence.date
    seurat.ETP@meta.data[['process.date']] = process.date
    
    #save basicQC features in the figure directory
    print("writing QC figure")
    png(paste0(directories$figure.basic.qc, "Basic_QC_", sampleID, ".png"))
    print(VlnPlot(seurat.ETP, features = c("nCount_RNA", "nFeature_RNA", "perc.mito", "nCount_ADT", "nFeature_ADT"), pt.size = 0))
    dev.off()
    
    #normalize data / find 50 variable features
    seurat.ETP = NormalizeData(seurat.ETP)
    seurat.ETP = FindVariableFeatures(seurat.ETP)
    
    #cell cycle regression
    cycle3 = fread('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/Signature_GeneList/regev_lab_cell_cycle_genes.txt', header = F)$V1
    s.genes = cycle3[1:43]
    g2m.genes = cycle3[44:97]
    seurat.ETP = CellCycleScoring(object = seurat.ETP, s.features = s.genes,
                                  g2m.features = g2m.genes, set.ident = TRUE)
    
    #heat shock regression
    heat_shock_gene = fread('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/Signature_GeneList/heat_shock_geneList.txt')
    heat_shock_gene = heat_shock_gene$`Approved symbol`
    seurat.ETP = AddModuleScore(seurat.ETP, features = list(heat_shock_gene),
                                name = 'HeatShock.Score')
    #ScaleData
    seurat.ETP = ScaleData(seurat.ETP, features = VariableFeatures(seurat.ETP), 
                           vars.to.regress = c('S.Score', 'G2M.Score', 'perc.mito',
                                               'HeatShock.Score1', 'nCount_RNA'))
    
    seurat.ETP = RunPCA(seurat.ETP, verbose = FALSE, npcs = 50)
    seurat.ETP = FindNeighbors(seurat.ETP, dims = 1:50)
    seurat.ETP = FindClusters(seurat.ETP, resolution = 1)
    
    ## remove doublets using DoubletFinder
    seurat.ETP = FindDoublets(seurat.ETP, exp_rate = 0.05)
    seurat.ETP = RunUMAP(seurat.ETP, dims = 1:50)
    p1 <- DimPlot(seurat.ETP, reduction = 'umap', 
                  group.by = 'Doublet_Singlet') + ggtitle('With Doublets')
    
    ## run basic seurat process after removing doublets
    sample.df$doublets[counter] = sum(seurat.ETP$Doublet_Singlet != 'Singlet')
    
    seurat.ETP = subset(seurat.ETP, Doublet_Singlet == 'Singlet')
    seurat.ETP = NormalizeData(seurat.ETP)
    seurat.ETP = FindVariableFeatures(seurat.ETP)
    seurat.ETP = ScaleData(seurat.ETP, features = VariableFeatures(seurat.ETP), 
                           vars.to.regress = c('S.Score', 'G2M.Score', 'perc.mito',
                                               'HeatShock.Score1', 'nCount_RNA'))
    seurat.ETP = RunPCA(seurat.ETP, verbose = FALSE, npcs = 50)
    seurat.ETP = RunUMAP(seurat.ETP, dims = 1:50)
    seurat.ETP = RunTSNE(seurat.ETP, dims = 1:50)
    seurat.ETP = FindNeighbors(seurat.ETP, reduction = 'pca', dims = 1:50)
    seurat.ETP = FindClusters(seurat.ETP, resolution = 1)
    
    p2 <- DimPlot(seurat.ETP, reduction = 'umap', label = T) + ggtitle('Doublets Removed')
    
    ggsave(CombinePlots(plots = list(p1, p2)), device = 'eps', width = 14, 
           height = 6, filename = paste0(directories$figure.doublet.plot, sampleID, '_umap_with_doublets.eps'))
    
    #normalize ADT
    seurat.ETP <- NormalizeData(seurat.ETP, assay = "ADT", normalization.method = "CLR")
    seurat.ETP <- ScaleData(seurat.ETP, assay = "ADT")
    
    sample.df$final.cells[counter] = ncol(seurat.ETP)
    
    #save seurat object----
    
    saveRDS(seurat.ETP, paste0(directories$seurat.object.directory, sampleID, '_CITESeq.rds'))

    print(paste0("saved seurat object to ", directories$seurat.object.directory, sampleID, '_CITESeq.rds'))
    
    print(paste("done for sample", sampleID, "moving onto next sample", counter+1, "out of", length(sample.df$sample.name)))
    
    
  }else{
    
    data <-  Read10X(data.dir = sample.df$matrix.location[counter])
    sampleID <- sample.df$sample.name[counter] %>% print()
    sample.group <- sample.df$group[counter] %>% print()
    sequence.date <- sample.df$sequence.date[counter]
    process.date <- sample.df$process.date[counter]
    
    print(paste("processing sample", sampleID, "which is", counter, "out of", length(sample.df$sample.name), "samples"))
    
    
    colnames(data) = paste0(sampleID, '_', colnames(data))
    #colnames(data$`Antibody Capture`) = paste0(sampleID, '_', colnames(data$`Antibody Capture`))
    
    original.dim <- dim(data) %>% print()
    
    sample.df$starting.cells[counter] = ncol(data)
    sample.df$starting.genes[counter] = nrow(data)
    sample.df$hb.filter[counter] = sum(data['HBB', ] > 3)
    
    print("starting to filter data")
    
    #filter hemoglobin----
    if('HBB' %in% rownames(data)){
      hbb_exp = data['HBB', ]
      data = data[, hbb_exp < 3]
    }
    
    if(sum(hbb_exp>3)>0){
      filter.rbc <- table(hbb_exp <3)
      print(paste("filtered", filter.rbc[1], "rbcs out, kept", filter.rbc[2], "cells"))
    }else{
      filter.rbc <- c(0, table(hbb_exp <3))
      print(paste("filtered", filter.rbc[1], "rbcs out, kept", filter.rbc[2], "cells"))
    }
    #filter by count between 1500 and 40,000
    nCount = Matrix::colSums(data) #totalreads per cell
    
    sample.df$nCount_1500.filter[counter] = sum(nCount<1500)
    sample.df$nCount40000.filter[counter] = sum(nCount>40000)
    
    
    data = data[, nCount < 40000 & nCount > 1500]
    
    sample.df$median.nCount[counter] = median(Matrix::colSums(data))
    
    filter.nCount <- table(nCount >1500 & nCount <40000) 
    print(paste("filtered", filter.nCount[1], "cells out, kept", filter.nCount[2], "cells"))
    
    #filter by perc.mito > 0.1
    mito.features <- grep(pattern = "^MT-", 
                          x = rownames(x = data), value = TRUE)
    
    perc.mito = Matrix::colSums(data[mito.features, ])/Matrix::colSums(data)
    
    #plot all mito features
    #par(mfrow=c(3,5))
    #for (i in 1:length(mito.features)){
    #data[mito.features[i], ] %>% hist(main = mito.features[i], xlim = c(1,500))
    #}
    
    #hist(perc.mito, main = sampleID) 
    
    #quantile(perc.mito)
    
    filter.mito <- table(perc.mito <= 0.1) %>% print()
    
    sample.df$mito.filter[counter] = sum(perc.mito >= 0.1)
    
    data = data[, perc.mito <= 0.1]
    perc.mito = perc.mito[perc.mito <= 0.1]
    
    #Subset Antibody Capture cells
    #data$`Antibody Capture` = data$`Antibody Capture`[,colnames(data$`Antibody Capture`) %in% colnames(data)]
    
    #filter by ncount_adt < 2000
    #nCount_ADT <- Matrix::colSums(data$`Antibody Capture`)
    #data$`Antibody Capture` = data$`Antibody Capture`[, nCount_ADT < 2000]
    
    filter.adt <-c(NA,NA)
    sample.df$filter.adt.2000.hypothetical = NA
    #Subset Gene by Antibody counts
    
    #data = data[,colnames(data) %in% colnames(data$`Antibody Capture`)]
    
    #final filter stats
    filter.stats <- c(original.dim, filter.rbc, filter.nCount, filter.mito, filter.adt) %>% 
      matrix(ncol = 2, byrow = T, dimnames = list(c(paste(sampleID,"original.dim"), "filter.rbc", "filter.nCount", "filter.mito", "filter.adt(hypothetical)"), c("filtered", "kept"))) %>% print()
    
    print("writing filter statistics")
    write.csv(filter.stats, file = paste0(directories$filter.stats.directory, sampleID, "_filterstats"))
    
    
    
    #Create Seurat Object
    seurat.ETP = CreateSeuratObject(counts = data) #first gene exp
    #seurat.ETP[['ADT']] = CreateAssayObject(counts = data$`Antibody Capture`) #then ADT
    
    cnames = colnames(seurat.ETP) #cell names
    tmp.mito = data.table('perc' = perc.mito, 'cname' = names(perc.mito)) #datatable, perc mito and cells as column names
    setkey(tmp.mito, cname) #fixes key
    seurat.ETP@meta.data[['perc.mito']] = tmp.mito[cnames]$perc #set seurat object metadata for perc.mito
    
    #create metadata
    print("writing sample.name, sample.group, sequence.date, process.date, to sample metadata")
    seurat.ETP@meta.data[['sample.name']] = sampleID
    seurat.ETP@meta.data[['sample.group']] = sample.group
    seurat.ETP@meta.data[['sequence.date']] = sequence.date
    seurat.ETP@meta.data[['process.date']] = process.date
    
    #save basicQC features in the figure directory
    print("writing QC figure")
    png(paste0(directories$figure.basic.qc, "Basic_QC_", sampleID, ".png"))
    print(VlnPlot(seurat.ETP, features = c("nCount_RNA", "nFeature_RNA", "perc.mito"), pt.size = 0))
    dev.off()
    
    #normalize data / find 50 variable features
    seurat.ETP = NormalizeData(seurat.ETP)
    seurat.ETP = FindVariableFeatures(seurat.ETP)
    
    #cell cycle regression
    cycle3 = fread('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/Signature_GeneList/regev_lab_cell_cycle_genes.txt', header = F)$V1
    s.genes = cycle3[1:43]
    g2m.genes = cycle3[44:97]
    seurat.ETP = CellCycleScoring(object = seurat.ETP, s.features = s.genes,
                                  g2m.features = g2m.genes, set.ident = TRUE)
    
    #heat shock regression
    heat_shock_gene = fread('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/Signature_GeneList/heat_shock_geneList.txt')
    heat_shock_gene = heat_shock_gene$`Approved symbol`
    seurat.ETP = AddModuleScore(seurat.ETP, features = list(heat_shock_gene),
                                name = 'HeatShock.Score')
    #ScaleData
    seurat.ETP = ScaleData(seurat.ETP, features = VariableFeatures(seurat.ETP), 
                           vars.to.regress = c('S.Score', 'G2M.Score', 'perc.mito',
                                               'HeatShock.Score1', 'nCount_RNA'))
    
    seurat.ETP = RunPCA(seurat.ETP, verbose = FALSE, npcs = 50)
    seurat.ETP = FindNeighbors(seurat.ETP, dims = 1:50)
    seurat.ETP = FindClusters(seurat.ETP, resolution = 1)
    
    ## remove doublets using DoubletFinder
    seurat.ETP = FindDoublets(seurat.ETP, exp_rate = 0.05)
    seurat.ETP = RunUMAP(seurat.ETP, dims = 1:50)
    p1 <- DimPlot(seurat.ETP, reduction = 'umap', 
                  group.by = 'Doublet_Singlet') + ggtitle('With Doublets')
    
    ## run basic seurat process after removing doublets
    sample.df$doublets[counter] = sum(seurat.ETP$Doublet_Singlet != 'Singlet')
    
    seurat.ETP = subset(seurat.ETP, Doublet_Singlet == 'Singlet')
    seurat.ETP = NormalizeData(seurat.ETP)
    seurat.ETP = FindVariableFeatures(seurat.ETP)
    seurat.ETP = ScaleData(seurat.ETP, features = VariableFeatures(seurat.ETP), 
                           vars.to.regress = c('S.Score', 'G2M.Score', 'perc.mito',
                                               'HeatShock.Score1', 'nCount_RNA'))
    seurat.ETP = RunPCA(seurat.ETP, verbose = FALSE, npcs = 50)
    seurat.ETP = RunUMAP(seurat.ETP, dims = 1:50)
    seurat.ETP = RunTSNE(seurat.ETP, dims = 1:50)
    seurat.ETP = FindNeighbors(seurat.ETP, reduction = 'pca', dims = 1:50)
    seurat.ETP = FindClusters(seurat.ETP, resolution = 1)
    
    p2 <- DimPlot(seurat.ETP, reduction = 'umap', label = T) + ggtitle('Doublets Removed')
    
    ggsave(CombinePlots(plots = list(p1, p2)), device = 'eps', width = 14, 
           height = 6, filename = paste0(directories$figure.doublet.plot, sampleID, '_umap_with_doublets.eps'))
    
    #normalize ADT
    #seurat.ETP <- NormalizeData(seurat.ETP, assay = "ADT", normalization.method = "CLR")
    #seurat.ETP <- ScaleData(seurat.ETP, assay = "ADT")
    
    sample.df$final.cells[counter] = ncol(seurat.ETP)
    
    #save seurat object----
    
    saveRDS(seurat.ETP, paste0(directories$seurat.object.directory, sampleID, '_RNASeq.rds'))
    
    print(paste0("saved seurat object to ", directories$seurat.object.directory, sampleID, '_CITESeq.rds'))
    
    print(paste("done for sample", sampleID, "moving onto next sample", counter+1, "out of", length(sample.df$sample.name)))
    
    
  }
  
}

#system(paste('rm -dr', directories$directories.list))
system(paste('mkdir -p', directories$directories.list))
saveRDS(directories, paste0(directories$directories.list, "directories.list_", number_samples, "samples_RNAseq_",Sys.Date(), '.rds'))
list.files(directories$directories.list)



