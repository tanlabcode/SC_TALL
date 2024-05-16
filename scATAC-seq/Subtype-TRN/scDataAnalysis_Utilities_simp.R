
## list of functions ####
library(data.table)
library(Matrix)
library(compiler)
library(magrittr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(heatmaply)
library(clusterProfiler)
library(reldist)
library(DoubletFinder)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)

library(mclust)
library(patchwork)
## **** note: for cluster validation, check:  ***** ##
## **** http://www.sthda.com/english/wiki/wiki.php?id_contents=7952  ***##


TF_IDF <- function (data, verbose = T) 
{
  if (class(x = data) == "data.frame") {
    data <- as.matrix(x = data)
  }
  if (class(x = data) != "dgCMatrix") {
    data <- as(object = data, Class = "dgCMatrix")
  }
  if (verbose) {
    message("Performing TF_IDF normalization")
  }
  npeaks <- colSums(x = data)
  tf <- t(x = t(x = data)/npeaks)
  idf <- ncol(x = data)/rowSums(x = data)
  norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
  norm.data[which(x = is.na(x = norm.data))] <- 0
  return(norm.data)
}


rebin_matrix2Bin <- function(mtx, resl = 100 * 1000){
  # mtx: matrix wiht rownames as chr-start-end and
  # colnames as cell names
  
  rnames = rownames(mtx)
  mtx_chr = sapply(rnames, function(x) unlist(strsplit(x, '-'))[1])
  chrs = unique(mtx_chr)
  starts = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[2]))
  ends = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[3]))
  
  rebin_mat = NULL
  for(chr0 in chrs){
    mtx0 = mtx[mtx_chr == chr0, ]
    mtx0 = as.data.table(mtx0)
    mtx0$start = starts[mtx_chr == chr0] 
    mtx0$end = ends[mtx_chr == chr0] 
    
    mtx0[, 'id' := ceiling((start+ end)/resl/2)]
    mtx0[, 'bin_id' := paste0(chr0, '-', id)]
    mtx0[, c('start', 'end', 'id') := NULL]
    rebin_mat = rbind(rebin_mat, mtx0)
  }
  
  rebin_mat = data.table(rebin_mat)
  setkey(rebin_mat, bin_id)
  
  new_mat = rebin_mat[, lapply(.SD, sum), by = bin_id]
  new_mat = new_mat[complete.cases(new_mat)]
  
  feature.names = new_mat$bin_id
  new_mat[, 'bin_id' := NULL]
  
  
  new_mat = as.matrix(new_mat)
  new_mat = as(new_mat, "sparseMatrix")
  rownames(new_mat) = feature.names
  
  return(new_mat)
}
rebin_matrix2Bin = cmpfun(rebin_matrix2Bin)

## given a matrix list, cbind them using the union set of features
cBind_union_features <- function(mat_list){
  ff = rownames(mat_list[[1]])
  for(i in 2:length(mat_list)){
    ff = unique(union(ff, rownames(mat_list[[i]])))
  }
  ## make a mtx with full features
  mat_union = list()
  for(i in 1:length(mat_list)){
    mtx0 = mat_list[[i]]
    ff0 = setdiff(ff, rownames(mtx0))
    if(length(ff0) > 0 ) {
      tmp = as(matrix(0, length(ff0), ncol(mtx0)), "sparseMatrix")
      rownames(tmp) = ff0
    }
    tmp_mat = rbind(mtx0, tmp)
    mat_union[[i]] = tmp_mat[order(rownames(tmp_mat)), ]
  }
  return(do.call('cbind', mat_union))
}
cBind_union_features = cmpfun(cBind_union_features)

# filtering of atac matrix
filterMat <- function(atac.mtx, minFrac_in_cell = 0.01, min_depth = 1000,
                      max_depth = 100000){
  depth.cell = Matrix::colSums(atac.mtx)
  atac.mtx = atac.mtx[, depth.cell > min_depth & depth.cell < max_depth]
  frac.in.cell = Matrix::rowSums(atac.mtx > 0)
  atac.mtx = atac.mtx[frac.in.cell > minFrac_in_cell, ]
  return(atac.mtx)
}


# assign gene to nearest peak and mark a gene if its tss within the peak
assignGene2Peak <- function(mtx, gene_ann, trans_dist = 1e+05){
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  peaks = tidyr::separate(data.table(x=rownames(mtx)),
                          col = x,
                          into = c('chr', 'start', 'end'))
  
  
  peaks$peak_name = rownames(mtx)
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  
  
  chrs = unique(peaks$chr)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    genes0 = gene_ann[chr == chr0]
    peaks0[, 'id' := which.min(abs(genes0$tss - start/2 - end/2)), by = peak_name]
    peaks0[, 'gene_name' := genes0[id, gene_name]]
    peaks0[, 'dist0' := min(abs(genes0$tss - start/2 -end/2)), by = peak_name]
    peaks0[, 'gene_name' := ifelse(dist0 > trans_dist, '', gene_name)]
    
    peaks0$tss_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= (peaks0$end[i] + 1000) & tss >= (peaks0$start[i] - 1000)] 
      if(nrow(tss0) > 0 ) {
        if(peaks0$gene_name[i] %in% tss0$gene_name) peaks0$gene_name[i] <- ''
        peaks0$tss_name[i] = paste(paste0(unique(tss0$gene_name), '-Tss'), 
                                                collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  peaks_ann[, 'id':= NULL]
  peaks_ann[, 'dist0':= NULL]
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(tss_name) & nchar(tss_name) > 1, 
                                         paste0(peak_new_name, ',', tss_name), peak_new_name)]
  setkey(peaks_ann, peak_name)
  
  
  
  rownames(mtx) = peaks_ann[rownames(mtx)]$peak_new_name
  
  return(mtx)
  
  
}


# assign gene to nearest peak and mark a gene if its tss within the peak
# input peak_coords with chr-start-end, format
assignGene2Peak_coords <- function(peak_coords, gene_ann, trans_dist = 1e+05){
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  peaks = tidyr::separate(data.table(x = peak_coords),
                          col = x,
                          into = c('chr', 'start', 'end'))
  
  
  peaks$peak_name = peak_coords
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  geneTssInPeak <- function(tss_ids, genes0){
    if(!is.na(tss))
      rr = genes0[tss <= end & tss >= start]$gene_name
    return(paste(rr, collapse = ','))
  }
  
  chrs = unique(peaks$chr)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    genes0 = gene_ann[chr == chr0]
    peaks0[, 'id' := which.min(abs(genes0$tss - start/2 - end/2)), by = 'peak_name']
    peaks0[, 'gene_name' := genes0[id, gene_name]]
    peaks0[, 'dist0' := min(abs(genes0$tss - start/2 -end/2)), by = peak_name]
    peaks0[, 'gene_name' := ifelse(dist0 > trans_dist, '', gene_name)]
    
    peaks0$tss_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= (peaks0$end[i] + 1000) & tss >= (peaks0$start[i] - 1000)]
      if(nrow(tss0) > 0 ) {
        if(peaks0$gene_name[i] %in% tss0$gene_name) peaks0$gene_name[i] <- ''
        peaks0$tss_name[i] = paste(paste0(unique(tss0$gene_name), '-Tss'), 
                                   collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  peaks_ann[, 'id':= NULL]
  peaks_ann[, 'dist0':= NULL]
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(tss_name) & nchar(tss_name) > 1, 
                                        paste0(peak_new_name, ',', tss_name), peak_new_name)]
  setkey(peaks_ann, peak_name)
  
  return(peaks_ann[peak_coords, ]$peak_new_name)
  
}


regress_on_pca <- function(seurat.obj, reg.var = 'nCount_ATAC'){
  
  pcs = seurat.obj@reductions$pca@cell.embeddings
  pcs.reg = pcs
  for(i in 1:length(reg.var)){
    
    reg.var0 = seurat.obj[[reg.var[i]]][[1]]
    pcs.reg = apply(pcs.reg, 2, function(x) lm(x ~ reg.var0)$residual )
    
  }
   colnames(pcs.reg) = colnames(pcs)
  seurat.obj@reductions$pca@cell.embeddings = pcs.reg
  return(seurat.obj)
}

# do normalization using log, tf-idf, or none, regress out confounds on pca or not 
doBasicSeurat_atac <- function(mtx, npc = 50, top.variable = 0.2, 
                               norm_by = c('log', 'tf-idf', 'none'),
                              doScale = T, doCenter = T, assay = 'ATAC',
                              reg.var = 'nCount_ATAC', regressOnPca = T,
                              project = 'scATAC'){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-')
 
  if(norm_by == 'log') seurat.obj@assays[[assay]]@data <- log1p(mtx) / log(2)
  if(norm_by == 'tf-idf') seurat.obj@assays[[assay]]@data <- TF.IDF(mtx, verbose = F)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = floor(nrow(mtx) * top.variable))
  
  if(regressOnPca){
    reg.var0 = NULL
  }else{
    reg.var0 = reg.var
  }
  
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = VariableFeatures(seurat.obj),
                          vars.to.regress = reg.var0, do.scale = doScale,
                          do.center = doCenter)
  
  
  #seurat.obj <- RunPCA(object = seurat.obj,
  #                     features = VariableFeatures(object = seurat.obj),
  #                     verbose = FALSE, seed.use = 10, npc = npc)
  seurat.obj <- RunPCA(object = seurat.obj,
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, npc = npc)
  if(length(reg.var) > 0 & regressOnPca) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  return(seurat.obj)
}

doBasicSeurat_atac = cmpfun(doBasicSeurat_atac)


# do normalization using log, tf-idf, or none, regress out confounds on pca or not 
doBasicSeurat_atac_updated <- function(mtx, npc = 30, top.variable = 5000, 
                               norm_by = c('log', 'tf-idf', 'none'),
                               doScale = T, doCenter = T, assay = 'ATAC',
                               reg.var = 'nCount_ATAC', regressOnPca = T,
                               project = 'scATAC', vap.min.frac = 0,
                               meta.data = NULL){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-',
                                  meta.data = meta.data)
  
  if(norm_by == 'log') seurat.obj@assays[[assay]]@data <- log1p(mtx) / log(2)
  if(norm_by == 'tf-idf') seurat.obj@assays[[assay]]@data <- TF.IDF(mtx, verbose = F)
  
  nvap = ifelse(top.variable > 1, top.variable, floor(top.variable * ncol(mtx)))
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = nvap)
  vaps = VariableFeatures(seurat.obj)
  peak.frac = rowMeans(mtx > 0)
  excludePks.fromVAP = names(which(peak.frac < vap.min.frac))
  vaps = setdiff(vaps, excludePks.fromVAP)
  
  if(length(vaps) < 10) stop('Top few VAPs left!')
  ## redo normalization using vap
  if(norm_by == 'tf-idf'){
    mtx.norm = TF.IDF(mtx[vaps, ])
    tmp <- mtx[setdiff(rownames(mtx), vaps), ]
    data0 <- rbind(mtx.norm, tmp)
    seurat.obj[[assay]]@data = data0[rownames(mtx), ]
    rm(data0, tmp, mtx.norm)
  }
  
  
  if(regressOnPca){
    reg.var0 = NULL
  }else{
    reg.var0 = reg.var
  }
  VariableFeatures(seurat.obj) <- vaps
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = vaps,
                          vars.to.regress = reg.var0, do.scale = doScale,
                          do.center = doCenter)
  
   seurat.obj <- RunPCA(object = seurat.obj,
                       features = vaps,
                       verbose = FALSE, npc = npc)
  if(length(reg.var) > 0 & regressOnPca) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  return(seurat.obj)
}

doBasicSeurat_atac_updated = cmpfun(doBasicSeurat_atac_updated)

#could be normalized by log, tf-idf or none
runSeurat_Atac <- function(mtx, npc = 50, top_variable_features = 5000, 
                           doScale = T, doCenter = T, assay = 'ATAC',
                           reg.var = NULL, norm_by = 'log', project = 'scATAC'){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-', min.cells = 1,
                                  min.features = 1)
  if(norm_by == 'log') seurat.obj[[assay]]@data <- log1p(seurat.obj[[assay]]@counts)/log(2)
  if(norm_by == 'tf-idf') seurat.obj[[assay]]@data <- TF_IDF(seurat.obj[[assay]]@counts)
  nvap = ifelse(top_variable_features > 1, top_variable_features, 
                floor(nrow(mtx) * top_variable_features))
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = nvap)
  
  ## remove variable features only accessible in less than 1% of cells
  mtx = seurat.obj[[assay]]@counts
  rs = Matrix::rowMeans(mtx > 0)
  rare.features = names(which(rs < 0.01))
  vaps = VariableFeatures(seurat.obj)
  vaps = setdiff(vaps, rare.features)
  niter = 0
  while(length(vaps) < 500 & nvap > 500){
    niter = niter + 1
    nvap = nvap + 2000
    seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                       selection.method = 'vst',
                                       nfeatures = min(nvap, nrow(seurat.obj)))
    vaps = VariableFeatures(seurat.obj)
    vaps = setdiff(vaps, rare.features)
    if(niter > 5) break
  }
  VariableFeatures(seurat.obj) <- vaps
  
  ## redo normalization using vap if norm by tf-idf
  if(norm_by == 'tf-idf'){
    mtx.norm = TF_IDF(mtx[vaps, ])
    tmp <- mtx[setdiff(rownames(mtx), vaps), ]
    data0 <- rbind(mtx.norm, tmp)
    seurat.obj[[assay]]@data = data0[rownames(mtx), ]
    rm(data0, tmp, mtx.norm)
  }
  
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = VariableFeatures(seurat.obj),
                          vars.to.regress = NULL, do.scale = doScale,
                          do.center = doCenter)
  
  
  seurat.obj <- RunPCA(object = seurat.obj,
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  if(length(reg.var) > 0 ) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  
  return(seurat.obj)
}
runSeurat_Atac = cmpfun(runSeurat_Atac)

# do normalization, pca using Seurat
doBasicSeurat_RNA <- function(mtx, npc = 50, top.variable = 0.2, pmito.upper = 0.2,
                           doScale = T, doCenter = T, min_nCount = 1500,
                           max_nCount = 15000, reg.var = 'nCount_RNA', sct=F,
                           mdata = NULL){

 ## top.variabl -- use top most variable features

 # filter cells with high percentage of mitocondria genes

  nCount = Matrix::colSums(mtx)
  mtx = mtx[, nCount < max_nCount & nCount > min_nCount]
  
  mito.features <- grep(pattern = "^MT-", 
                      x = rownames(x = mtx), value = TRUE)

  perc.mito = Matrix::colSums(mtx[mito.features, ])/Matrix::colSums(mtx)

  mtx = mtx[, perc.mito <= pmito.upper]
  perc.mito = perc.mito[perc.mito <= pmito.upper]


 # create seurat object
  if(!is.null(mdata)) mdata = mdata[colnames(mtx), ]
  seurat.obj = CreateSeuratObject(mtx, project = 'scRNA', assay = 'RNA',
                                  names.delim = '-', min.cells = 10, min.features = 100,
                                  meta.data = mdata)
  

  # add perc.mito to seurat objects
  cnames = colnames(seurat.obj)
  tmp.mito = data.table('perc' = perc.mito, 'cname' = names(perc.mito))
  setkey(tmp.mito, cname)
  seurat.obj@meta.data[['perc.mito']] = tmp.mito[cnames]$perc
  
  
  seurat.obj <- subset(x = seurat.obj, subset = (nFeature_RNA < 10000))
  vegs = ifelse(top.variable > 1, top.variable, floor(nrow(mtx) * top.variable))
  
  if(!sct){
    seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
                                scale.factor = 1e4)
    seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                       selection.method = 'vst', 
                                       nfeatures = vegs)
    seurat.obj <- ScaleData(object = seurat.obj, 
                            features = VariableFeatures(seurat.obj), 
                            vars.to.regress = reg.var, 
                            do.scale = doScale, do.center = doCenter)
    
    
  }else{
    seurat.obj <- SCTransform(seurat.obj, vars.to.regress = reg.var, verbose = F,
                              variable.features.n = floor(nrow(mtx) * top.variable))
  }
  
  #seurat.obj <- RunPCA(object = seurat.obj, 
  #                     features = VariableFeatures(object = seurat.obj),
  #                     verbose = FALSE, seed.use = 10, npc = npc)
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE,  npc = npc)
  
  return(seurat.obj)
}
doBasicSeurat_RNA = cmpfun(doBasicSeurat_RNA)



# Find doublets
FindDoublets <- function(seurat.rna, PCs = 1:50, 
                         exp_rate = 0.02, sct = FALSE){
  # sct--do SCTransform or not
  
  ## pK identification
  sweep.res.list <- paramSweep_v3(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25,
                                 pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = sct)
  
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25, 
                                 pK = 0.09, nExp = nExp_poi.adj,
                                 reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                                 sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  #seurat.rna = subset(seurat.rna, Doublet_Singlet == 'Singlet')
  return(seurat.rna)
}


doSeurat_rmDoublets <- function(sampleID, exp_rate = 0.05, pmito.upper = 0.15,
                                min_nCount = 1500, max_nCount = 50000,
                                clusterOn = 'pca', npc = 50,
                                resolution = 0.5){
  
  dir0 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/CellRangerResults/June19_2019/'
  
  sampleName = paste0('MLL_', sampleID, '_scRNA')
  mtx = Read10X(paste0(dir0, sampleName, '/', sampleName, '/outs/filtered_feature_bc_matrix'))
  colnames(mtx) = paste0(sampleID, '_', colnames(mtx))
  
  # remove red blood cells
  if('HBB' %in% rownames(mtx)){
    hbb_exp = mtx['HBB', ]
    mtx = mtx[, hbb_exp < 3]
  }
  
  
  ##try remove doublets instead of manully filter cells with larger UMI
  seurat.obj = doBasicSeurat_RNA(mtx, min_nCount = min_nCount, max_nCount = max_nCount,
                                 npc = npc, pmito.upper = pmito.upper)
  seurat.obj = FindNeighbors(seurat.obj, reduction = clusterOn, dims = 1:npc)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  seurat.obj = FindDoublets(seurat.obj, exp_rate = exp_rate)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  p1 <- DimPlot(seurat.obj, reduction = 'umap', 
                group.by = 'Doublet_Singlet') + ggtitle('With Doublets')
  
  ## remove doublets and do PCA and clustering again
  seurat.obj = subset(seurat.obj, Doublet_Singlet == 'Singlet')
  seurat.obj = doBasicSeurat_RNA(seurat.obj@assays$RNA@counts, 
                                 min_nCount = min_nCount, max_nCount = max_nCount)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  seurat.obj = RunTSNE(seurat.obj, dims = 1:npc)
  seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:npc)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  
  p2 <- DimPlot(seurat.obj, reduction = 'umap', label = T) + ggtitle('Doublets Removed')
  
  
  
  
  ## load to viscello
  inputs = prepInput4Cello(mtx = seurat.obj@assays$RNA@counts, 
                           seurat.obj = seurat.obj,
                           cello.name = sampleName,
                           assay = 'RNA', 
                           extraDR = T, cluster_excluded = NULL)
  
  # save inputs for future use
  celloInputDir = paste0('Input4VisCello/', sampleName)
  system(paste0('mkdir -p ', celloInputDir))
  saveRDS(inputs$eset, file = paste0(celloInputDir, '/eset.rds'))
  saveRDS(inputs$clist, file = paste0(celloInputDir, '/clist.rds'))
  
  ggsave(CombinePlots(plots = list(p1, p2)), device = 'eps', width = 14, 
         height = 6, filename = paste0('Figures/scRNA/MLL_', 
                                       sampleID, '/', sampleName, '_umap_with_doublets.eps'))
  
  
  saveRDS(seurat.obj, paste0('Seurat_Objects/scRNA/seurat_', sampleName, '_doubletRemoved.rds'))
  
  return(seurat.obj)
}

##bcPrefix can be set as sampleName or ID
doSeurat_rmDoublets_dir <- function(filtered_mtx_dir, exp_rate = 0.05, pmito.upper = 0.15,
                                min_nCount = 1500, max_nCount = 50000,
                                clusterOn = 'pca', npc = 50,
                                resolution = 0.5, bcPrefix = 'pbmc',
                                savePlot = T){
  
  
  mtx = Read10X(filtered_mtx_dir)
  colnames(mtx) = paste0(bcPrefix, '_', colnames(mtx))
  
  # remove red blood cells
  if('HBB' %in% rownames(mtx)){
    hbb_exp = mtx['HBB', ]
    mtx = mtx[, hbb_exp < 3]
  }
  
  
  ##try remove doublets instead of manully filter cells with larger UMI
  seurat.obj = doBasicSeurat_RNA(mtx, min_nCount = min_nCount, max_nCount = max_nCount,
                                 npc = npc, pmito.upper = pmito.upper)
  seurat.obj = FindNeighbors(seurat.obj, reduction = clusterOn, dims = 1:npc)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  seurat.obj = FindDoublets(seurat.obj, exp_rate = exp_rate)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  p1 <- DimPlot(seurat.obj, reduction = 'umap', 
                group.by = 'Doublet_Singlet') + ggtitle('With Doublets')
  
  ## remove doublets and do PCA and clustering again
  seurat.obj = subset(seurat.obj, Doublet_Singlet == 'Singlet')
  seurat.obj = doBasicSeurat_RNA(seurat.obj@assays$RNA@counts, 
                                 min_nCount = min_nCount, max_nCount = max_nCount)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:npc)
  seurat.obj = RunTSNE(seurat.obj, dims = 1:npc)
  seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:npc)
  seurat.obj = FindClusters(seurat.obj, resolution = resolution)
  p2 <- DimPlot(seurat.obj, reduction = 'umap', label = T) + 
    ggtitle('Doublets Removed')
  
  if(savePlot) ggsave(CombinePlots(plots = list(p1, p2)), device = 'eps', width = 14, 
         height = 6, filename = paste0('Figures/scRNA/', bcPrefix, '_umap_with_doublets.eps'))
  
  seurat.obj$sample = bcPrefix
  return(seurat.obj)
}

# Find doublets
FindDoublets_Atac <- function(seurat.atac, PCs = 1:50, 
                         exp_rate = 0.02, sct = FALSE){
  # sct--do SCTransform or not
  
  if(!sct){
    seurat.rna = CreateSeuratObject(seurat.atac@assays$ATAC@counts)
    seurat.rna = NormalizeData(seurat.rna)
    seurat.rna = FindVariableFeatures(seurat.rna)
    VariableFeatures(seurat.rna) <- VariableFeatures(seurat.atac)
    seurat.rna = ScaleData(seurat.rna)
    seurat.rna = RunPCA(seurat.rna, npcs = max(PCs), verbose = F)
    seurat.rna@reductions$pca@cell.embeddings <- seurat.atac@reductions$pca@cell.embeddings
    seurat.rna@reductions$pca@feature.loadings <- seurat.atac@reductions$pca@feature.loadings
    seurat.rna$seurat_clusters = seurat.atac$seurat_clusters
  }else{
    seurat.rna = seurat.atac
    seurat.rna@assays$RNA <- seurat.atac@assays$ATAC
  }
  
  ## pK identification
  sweep.res.list <- paramSweep_v3(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25,
                                 pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = sct)
  
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25, 
                                 pK = 0.09, nExp = nExp_poi.adj,
                                 reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                                 sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  #seurat.rna = subset(seurat.rna, Doublet_Singlet == 'Singlet')
  seurat.atac$Doublet_Singlet = seurat.rna$Doublet_Singlet
  return(seurat.atac)
}




# basic seurat plot
basicSeuratPlot <- function(org.seurat.obj){
  p1 <- VizDimLoadings(object = org.seurat.obj, dims = 1:2, nfeatures = 20)
  
  p2 <- DimHeatmap(object = org.seurat.obj, dims = 1:6, cells = 100,   balanced = TRUE)
  p3 <- ElbowPlot(object = org.seurat.obj, ndims = 50)
  mfeatures = names(org.seurat.obj@meta.data)
  pfeatures = mfeatures[grepl(mfeatures, pattern = 'nCount_|mito')]
  p4 <- VlnPlot(org.seurat.obj, pfeatures)
  return(list(p1, p2, p3, p4))
}


## map gene to atac peak
gene2peak <- function(gene_set, peaks, gene_ann){
  # should include tss information in gene_list
  gene_list = gene_ann[gene_name %in% gene_set, ]
  chrs = unique(gene_list$chr)
  gene_new = NULL
  peaks[, 'midP' := start/2 + end/2]
  for(chr0 in chrs){
    gene0 = gene_list[chr == chr0, ]
    peaks0 = peaks[chr == chr0]
    gene0[, 'peak_id' := which(tss >= peaks0$start & tss <= peaks0$end), by = gene_id]
    gene0[, 'peak_id' := ifelse(is.na(peak_id), which.min(abs(tss - peaks0$midP)), peak_id), by = gene_name]
    gene0[, 'peak' := peaks0[peak_id]$pos]
    
    gene_new = rbind(gene_new, gene0)
  }
  return(gene_new)
}


read10X_ATAC <- function(dirt){
  mtx_path <- paste0(dirt, "matrix.mtx")
  feature_path <- paste0(dirt, "peaks.bed")
  barcode_path <- paste0(dirt, "barcodes.tsv")
  
  
  features <-readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature, sep = '-')
  barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
  
  mtx <-  Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$feature)%>%
    magrittr::set_colnames(barcodes$barcode) 
  
  return(mtx)
}

read_mtx_scATACpro <- function(mtx_path){
  #mtx_path <- paste0(dirt, "matrix.mtx")
  mtx.dir = dirname(mtx_path)
  feature_path <- paste0(mtx.dir, "/features.txt")
  barcode_path <- paste0(mtx.dir, "/barcodes.txt")
  
  
  features <- fread(feature_path, header = F)
  barcodes <- fread(barcode_path, header = F)
  
  mtx <-  Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$V1)%>%
    magrittr::set_colnames(barcodes$V1)
  
  return(mtx)
}


normalize_gene_activities.corrected <- function (activity_matrices, cell_num_genes) 
{
  if (!is.list(activity_matrices)) {
    scores <- activity_matrices
    normalization_df <- data.frame(cell = colnames(activity_matrices), 
                                   cell_group = 1)
  }else {
    scores <- do.call(cbind, activity_matrices)
    normalization_df <- do.call(rbind, lapply(seq_along(activity_matrices), 
                                              function(x) {
                                                data.frame(cell = colnames(activity_matrices[[x]]), 
                                                           cell_group = rep(x, ncol(activity_matrices[[x]])))
                                              }))
  }
  scores <- scores[Matrix::rowSums(scores) != 0, Matrix::colSums(scores) != 
                     0]
  
  ## correct by adding following lines
  cell_num_genes = cell_num_genes[normalization_df$cell %in% colnames(scores)]
  normalization_df = subset(normalization_df, cell %in% colnames(scores))
  ##
  
  normalization_df$cell_group <- factor(normalization_df$cell_group)
  normalization_df$total_activity <- Matrix::colSums(scores)
  normalization_df$total_sites <- cell_num_genes[as.character(normalization_df$cell)]
  if (!is.list(activity_matrices)) {
    activity_model <- stats::lm(log(total_activity) ~ log(total_sites), 
                                data = normalization_df)
  } else {
    activity_model <- stats::lm(log(total_activity) ~ log(total_sites) * 
                                  cell_group, data = normalization_df)
  }
  normalization_df$fitted_curve <- exp(as.vector(predict(activity_model, 
                                                         type = "response")))
  size_factors <- log(normalization_df$fitted_curve)/mean(log(normalization_df$fitted_curve))
  size_factors <- Matrix::Diagonal(x = 1/size_factors)
  row.names(size_factors) <- normalization_df$cell
  colnames(size_factors) <- row.names(size_factors)
  
  scores <- Matrix::t(size_factors %*% Matrix::t(scores))
  scores@x <- pmin(1e+09, exp(scores@x) - 1)
  sum_activity_scores <- Matrix::colSums(scores)
  scale_factors <- Matrix::Diagonal(x = 1/sum_activity_scores)
  row.names(scale_factors) <- normalization_df$cell
  colnames(scale_factors) <- row.names(scale_factors)
  scores <- Matrix::t(scale_factors %*% Matrix::t(scores))
  
  if (!is.list(activity_matrices)) {
    rn = row.names(activity_matrices)
    cn = colnames(activity_matrices)
    rn = rn[rn %in% row.names(scores)]
    cn = cn[cn %in% colnames(scores)]
    
    ret <- scores[rn, cn]
  }
  else {
    ret <- lapply(activity_matrices, function(x) {
      rn = row.names(x)
      cn = colnames(x)
      rn = rn[rn %in% row.names(x)]
      cn = cn[cn %in% colnames(x)]
      scores[rn, cn]
    })
  }
  return(ret)
}

# do cicero given a Seurat object, output gene activity score
doCicero_gascore <- function(seurat.obj, reduction = 'umap', chr_sizes,
                             gene_ann, npc = 30, coaccess_thr = 0.25,
                             min_frac_cell = 0.01){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019)
  library(cicero)
  mtx = GetAssayData(seurat.obj, slot = 'counts')
  mtx = 1*(mtx > 0)
  ## remove peaks that only appear in less than 0.5% of cells
  rs = Matrix::rowMeans(mtx > 0)
  mtx = mtx[rs > min_frac_cell, ]
  
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) unlist(strsplit(x, ','))[1])
  
  new.rnames = sapply(new.rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  
  
  #dt = reshape2::melt(as.matrix(mtx), value.name = 'Count')
  #dt = dt[dt$Count > 0, ]
  dt = mefa4::Melt(mtx)
  rm(mtx)
  names(dt) = c('Peak', 'Cell', 'Count')
  
  
  dt$Cell = as.character(dt$Cell)
  dt$Peak = as.character(dt$Peak)
  
  input_cds <- make_atac_cds(dt, binarize = T)
  rm(dt)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  if(reduction == 'tsne') {
    if(is.null(seurat.obj@reductions$tsne))
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, check_duplicates = F)
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    if(is.null(seurat.obj@reductions$umap))
      seurat.object <- RunUMAP(seurat.object, dims = 1:npc)
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  #make the cell id consistet
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords)
  
  ## get connections
  
  conns <- run_cicero(cicero_cds, chr_sizes)
  
  ## get cicero gene activity score
  names(gene_ann)[4] <- "gene"
  
  input_cds <- annotate_cds_by_site(input_cds, gene_ann)
  
  # generate unnormalized gene activity matrix
  unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
  
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  
  # normalize
  cicero_gene_activities <- normalize_gene_activities.corrected(unnorm_ga, num_genes)
  
  # if you had two datasets to normalize, you would pass both:
  # num_genes should then include all cells from both sets
  #unnorm_ga2 <- unnorm_ga
  #cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), num_genes)
  conns = data.table(conns)
  conns = conns[coaccess > coaccess_thr, ]
  res = list('conns' = conns, 'ga_score' = cicero_gene_activities)
  return(res)  
}


# do cicero given a Seurat object, just return the connection 
doCicero_conn <- function(seurat.obj, reduction = 'tsne', 
                          chr_sizes, npc = 30, coaccess_thr = 0.25,
                          min_frac_cell = 0.01){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019)
  library(cicero)
  mtx = GetAssayData(seurat.obj, slot = 'counts')
  rs = Matrix::rowMeans(mtx > 0)
  mtx = mtx[rs > min_frac_cell, ]
  
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) unlist(strsplit(x, ','))[1])
  new.rnames = sapply(new.rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  dt = mefa4::Melt(mtx)
  rm(mtx)
  
  names(dt) = c('Peak', 'Cell', 'Count')
  dt$Cell = as.character(dt$Cell)
  dt$Peak = as.character(dt$Peak)
  input_cds <- make_atac_cds(dt, binarize = T)
  rm(dt)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  if(reduction == 'tsne') {
    if(is.null(seurat.obj@reductions$tsne))
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, check_duplicates = F)
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    if(is.null(seurat.obj@reductions$umap))
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc)
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  #make the cell id consistet
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords)
  
  ## get connections
  
  conns <- run_cicero(cicero_cds, chr_sizes)
  conns = data.table(conns)
  conns = conns[coaccess > coaccess_thr, ]
  return(conns)
}


run_chromVAR <- function(mtx, genomeName = 'BSgenome.Hsapiens.UCSC.hg38',
                         ncore = 3){
  
  register(MulticoreParam(ncore))
  if(!require(genomeName, character.only = T)) BiocManager::install(genomeName)
  
  peaks = data.table('x' = rownames(mtx))
  peaks = tidyr::separate(peaks, col = 'x', into = c('chr', 'start', 'end'), sep = '-')
  peaks = GenomicRanges::makeGRangesFromDataFrame(peaks)
  
  frag.counts = SummarizedExperiment(assay = list(counts = mtx),
                                     rowRanges = peaks)
  frag.counts <- addGCBias(frag.counts, genome = genomeName)
  #motifs <- getJasparMotifs()
  library(chromVARmotifs) ## cisbp motif
 
  if(grepl(genomeName, pattern = 'hg')){
    motifs = human_pwms_v2
  }else{
    motifs = mouse_pwms_v2
  }
  
  motif_ix <- matchMotifs(motifs, frag.counts,
                          genome = genomeName)
  dev <- computeDeviations(object = frag.counts, 
                           annotations = motif_ix)
  bg <- getBackgroundPeaks(object = frag.counts)
  
  dev <- computeDeviations(object = frag.counts, annotations = motif_ix,
                           background_peaks = bg)
  
  #motif.zscore = dev@assays$data$z
  return(dev)
}



## get overlap number of nearby DAs of DE within a given dist
getOverlap_DEwithDA <- function(DEs, DAs, flank.dist = 500 * 1000){
  chrs = unique(DEs$chr)
  res.overlap = NULL
  for(chr0 in chrs){
    DE0 = DEs[chr == chr0, ]
    DA0 = DAs[chr == chr0, ]
    if(nrow(DA0) == 0) next
    tmp.res = lapply(DE0$midP, function(x) DA0[abs(midP - x) <= flank.dist])
    tmp.res = do.call('rbind', tmp.res)
    res.overlap = rbind(res.overlap, tmp.res)
  }
  res.overlap = res.overlap[!duplicated(res.overlap), ]
  return(res.overlap)
}

## link DEs with DAs 
## given a table of DE and a tables DAs (with cluster info), output a matrix records the number of 
## DAs for each groups DEs with a flank.dist 
## the DE_table should include columns <cluster> <gene> <p_val>
## and DA_table should include columns <cluster> <peak> <p_val>
## Note: peak should be like chr-start-end, p_val or logFC should be 
linkDEwithDA <- function(DE_table, DA_table, gene_ann, flank.dist = 500*1000, 
                      top.de = 0.2, top.da = 0.2, sort_by = 'p_val', enhs = NULL){
  de.cl = unique(DE_table$cluster)
  da.cl = unique(DA_table$cluster)
  len1 = length(de.cl)
  len2 = length(da.cl)
  noverlaps = noverlaps.a2e = matrix(0, len1, len2)
 
  # add coordinate to DE
  gene_ann = gene_ann[!duplicated(gene_name), ]
  setkey(gene_ann, gene_name)
  DE_table = DE_table[gene %in% gene_ann$gene_name, ]
  DE_table[, 'chr' := gene_ann[J(DE_table$gene)]$chr]
  DE_table[, 'midP' := gene_ann[J(DE_table$gene)]$tss]
  
  # get midpoint of DA
  DA_table[, 'chr':= unlist(strsplit(peak, '-'))[1], by = peak]
  DA_table[, 'start':= as.integer(unlist(strsplit(peak, '-'))[2]), by = peak]
  DA_table[, 'end':= as.integer(unlist(strsplit(peak, '-'))[3]), by = peak]
  DA_table[, 'midP' := floor(start/2 + end/2)]
  
  # only use peak overlapped with enhancers
  if(!is.null(enhs)){
    filter.da = NULL
    chrs = unique(DA_table$chr)
    for(chr0 in chrs){
      enh0 = enhs[chr == chr0, ]
      DA0 = DA_table[chr == chr0, ]
      DA0[, 'mdist' := min(abs(midP - enh0$midP)), by = midP]
      DA0[, 'len' := end - start]
      
      filter.da = rbind(filter.da, DA0[mdist < 1000 + len/2])
    }
    filter.da[, c('mdist', 'len') := NULL]
    DA_table = filter.da
    rm(filter.da)
  }
  
  # calculate overlaps
  for(i in 1:len1){
    DEs = DE_table[cluster == de.cl[i], ]
    setkeyv(DEs, sort_by)
    for(j in 1:len2){
      DAs = DA_table[cluster == da.cl[j], ]
      setkeyv(DAs, sort_by)
      topDE = floor(nrow(DEs) * top.de)
      topDA = floor(nrow(DAs) * top.da)
      noverlaps[i, j] = nrow(getOverlap_DEwithDA(DEs[1:topDE, ], DAs[1:topDA, ],
                                            flank.dist))
      noverlaps.a2e[i, j] = nrow(getOverlap_DEwithDA(DAs[1:topDA, ], DEs[1:topDE, ], 
                                                flank.dist))
    }
  }
  
  # transform to pvalue(hyper-geometric test), given overlap mat and
  # col and row sums (union)
  get_pv <- function(noverlaps, da.count, de.count, nn){
    #da.count = colSums(noverlaps)
    #de.count = rowSums(noverlaps)
    
    da.mat = matrix(rep(da.count, each = nrow(noverlaps)), 
                    nrow = nrow(noverlaps))
    de.mat = matrix(rep(de.count, each = ncol(noverlaps)), 
                    nrow = nrow(noverlaps), byrow = T)
    
    pvs = phyper(noverlaps, da.mat, nn - da.mat,
                 de.mat, lower.tail = F)
    # normalized overlaps
    noverlaps.normDA = noverlaps/(da.mat * de.mat) * nn
    #tt = t(percentize(t(noverlaps.normDA)))
    
    return(list(pvs, noverlaps.normDA))
  }
  
 
  col.pool = lapply(da.cl, function(x) getOverlap_DEwithDA(DE_table, DA_table[cluster == x, ],  flank.dist))
  row.pool = lapply(de.cl, function(x) getOverlap_DEwithDA(DE_table[cluster == x, ], DA_table,  flank.dist))
  all.over = do.call('rbind', col.pool)
  all.over = all.over[!duplicated(all.over), ] ## should remove duplicates
  nn = nrow(all.over) # total 'balls'
  
  res = get_pv(noverlaps, sapply(col.pool, function(x) nrow(x)), 
               sapply(row.pool, function(x) nrow(x)), nn)
  
  col.pool = lapply(da.cl, function(x) getOverlap_DEwithDA(DA_table[cluster == x, ], DE_table, flank.dist))
  row.pool = lapply(de.cl, function(x) getOverlap_DEwithDA(DA_table, DE_table[cluster == x, ],  flank.dist))
  all.over = do.call('rbind', col.pool)
  all.over = all.over[!duplicated(all.over), ]
  nn = nrow(all.over)
  res.a2e = get_pv(noverlaps.a2e, sapply(col.pool, function(x) nrow(x)), 
                   sapply(row.pool, function(x) nrow(x)), nn)
  
  return(list('pvs.e2a' = res[[1]], 'norm.e2a' = res[[2]], 
               'pvs.a2e' = res.a2e[[1]], 'norm.a2e' = res.a2e[[2]],
              'count.e2a' = noverlaps, 'count.a2e' = noverlaps.a2e,
              'DA_filtered' = DA_table))
  
}


# mtx: gene by cell matrix, or peak/bin by cell matrix
# seurat.obj: optional, if not provided, will construct one (which will clustering on pca1:50),
# if seurat.obj is provided, assume it did pca
# norm_mtx: normalized matrix, equals to log(mtx +1) if not provided
# extraDims: do extra dimension reduction on 10, 30, 100  etc
# subSamples: provided a subsample version
# vFeatures: given variable features; if NULL using default variable features
prepInput4Cello <- function(mtx, seurat.obj, norm_mtx = NULL, 
                            cello.name = 'scRNA', assay = 'RNA', 
                            resl4clust = 0.6, 
                            extraDims = c(10, 20, 30, 50, 80, 100),
                            subSamples = NULL, subCelloName = 'sub',
                            vars.to.regOnPca = NULL, downSample = NULL,
                            vFeatures = NULL, reduction = 'pca'){

  if(!is.null(vFeatures)) seurat.obj <- ScaleData(seurat.obj, features = vFeatures)
  if(is.null(vFeatures))  vFeatures = VariableFeatures(seurat.obj)  
  DefaultAssay(seurat.obj) = assay
  
  if(!is.null(seurat.obj@reductions$tsne)) my_tsne_proj <- seurat.obj@reductions$tsne@cell.embeddings
  if(!is.null(seurat.obj@reductions$umap)) my_umap_proj <- seurat.obj@reductions$umap@cell.embeddings 
  if(!is.null(seurat.obj@reductions$pca)) {
    my_pca_proj <- seurat.obj@reductions$pca@cell.embeddings
  
    ## change pca name to be consistent with cello
    pc.name = colnames(my_pca_proj)
    colnames(my_pca_proj) <- sapply(pc.name, function(x) gsub('_', '', x))
    ndefault = ncol(my_pca_proj)
  }
  if(reduction == 'harmony') ndefault = ncol(seurat.obj@reductions$harmony@cell.embeddings)
  #mtx = mtx[rownames(mtx) %in% rownames(seurat.obj), ]
  mtx = mtx[, colnames(mtx) %in% colnames(seurat.obj)]
  
  meta = seurat.obj@meta.data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  if(is.null(norm_mtx)) norm_mtx = log2(mtx + 1)
  
  ids = 1:ncol(seurat.obj)
  if(!is.null(downSample)){
    set.seed(2019)
    if(downSample <= ncol(seurat.obj)) 
      ids = sample(1:ncol(seurat.obj), downSample)
  }
  
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx[, ids], 
                                       norm_exprs = norm_mtx[, ids]),
              phenoData =  new("AnnotatedDataFrame", data = meta[ids, ]),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  # Creating a cello 
  Cello <- setClass("Cello",
                    slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                    )
  )
  
  cello1 <- new("Cello", name = cello.name, idx = 1:length(ids)) 
  
  cello1@proj <- list("PCA" = my_pca_proj[ids, ]) 
  if(!is.null(seurat.obj@reductions$tsne)) cello1@proj[[paste0('t-SEN', 'Orig')]] <- my_tsne_proj[ids, ]
  if(!is.null(seurat.obj@reductions$umap)) cello1@proj[[paste0('UMAP', 'Orig')]] <- my_umap_proj[ids, ]
  
  
  if(is.null(extraDims)) extraDims = ndefault
  
  for(dim0 in extraDims){
    
    if(reduction == 'pca') seurat.obj <- RunPCA(seurat.obj, npcs = dim0, verbose = F,
                         assay = assay, features = vFeatures)
    
    if(!is.null(vars.to.regOnPca)) seurat.obj = regress_on_pca(seurat.obj, vars.to.regOnPca)
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:dim0, check_duplicates = FALSE,
                          assay = assay, reduction = reduction)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:dim0, verbose = F,
                          assay = assay, reduction = reduction)
    my_tsne_proj0 <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj0 <- seurat.obj@reductions$umap@cell.embeddings
    
    
    
    cello1@proj[[paste0('t-SEN', dim0)]] <- my_tsne_proj0[ids, ]
    cello1@proj[[paste0('UMAP', dim0)]] <- my_umap_proj0[ids, ]
   
  }
 
  #cello1@proj = cello1@proj[order(names(cello1@proj))]

  clist <- list()
  clist[[cello.name]] <- cello1
 
  if(!is.null(subSamples)){
    ids = which(colnames(seurat.obj) %in% subSamples)
    cello2 <- new("Cello", name = cello.name, idx = ids) 
    
    
    cello2@proj <- list("PCA" = my_pca_proj[ids, ], 
                        't-SNE5' = my_tsne_proj[ids, ],
                        "UMAP5" = my_umap_proj[ids, ]) 
   
    cello.name.sub = paste0(cello.name, '_', subCelloName)
    clist[[cello.name.sub]] <- cello2
  }
  
  
  return(list('eset' = eset, 'clist' = clist))
}

## update different version of clustering, and other information
## suppose cells are the same in seurat.obj, and only 1 clist obj
## updateProjDims means update tsne and umap using pca updateProjDims
## add extra meta data for each cell
updateInput4Cello <- function(celloInputPath, seuratObjPath, 
                              assay = 'integrated', extraMeta = NULL,
                              addPcaDims = c(20),
                              clusterOn = 'pca', clusterOnDim = 30,
                              resolutions = c(0.2, 0.4, 0.6)){
  eset = readRDS(paste0(celloInputPath, '/eset.rds'))
  clist = readRDS(paste0(celloInputPath, '/clist.rds'))
  mdata = phenoData(eset)
  
  seurat.obj = readRDS(seuratObjPath)
  DefaultAssay(seurat.obj) <- assay
  
  
  
  ## check whether clusterOnDim exists
  if(!is.null(clusterOnDim)){
    if(ncol(seurat.obj@reductions$pca@cell.embeddings) < clusterOnDim){
      seurat.obj <- RunPCA(seurat.obj, npcs = clusterOnDim, verbose = F)
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:clusterOnDim, check_duplicates = FALSE,
                            assay = assay)
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:clusterOnDim, verbose = F,
                            assay = assay)
      my_tsne_proj <- seurat.obj@reductions$tsne@cell.embeddings
      my_umap_proj <- seurat.obj@reductions$umap@cell.embeddings 
      
      clist[[1]]@proj[[paste0('t-SNE', clusterOnDim)]] = my_tsne_proj
      clist[[1]]@proj[[paste0('UMAP', clusterOnDim)]] = my_umap_proj
      
    }
    
  }
  
  ## try different clustering
  if(!is.null(clusterOnDim)){
    seurat.obj <- FindNeighbors(seurat.obj, reduction = clusterOn, dims = 1:clusterOnDim)
    for(res0 in resolutions){
      seurat.obj <- FindClusters(seurat.obj, resolution = res0)
      mdata@data[[paste0(assay, '_snn_', clusterOn, clusterOnDim , 'res.', res0)]] <- 
        seurat.obj[[paste0(assay, '_snn_res.', res0)]][, 1]
    }
  }
  
  mdata@data$seurat_clusters = seurat.obj$seurat_clusters
  
  if(!is.null(extraMeta)){
    if(all(rownames(extraMeta) == rownames(mdata))) {
      shared.features = intersect(colnames(mdata@data), colnames(extraMeta))
      mdata@data[, shared.features] <- NULL
      mdata@data = cbind(mdata@data, extraMeta)
    }
  }
  
  phenoData(eset) <- mdata
  
  
  ## update clist
  for(dim0 in addPcaDims){
    if(ncol(seurat.obj@reductions$pca@cell.embeddings) < dim0) {
      seurat.obj <- RunPCA(seurat.obj, verbose = F, npcs = dim0)
    }
    
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:dim0, check_duplicates = FALSE,
                          assay = assay)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:dim0, verbose = F,
                          assay = assay)
    my_tsne_proj <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj <- seurat.obj@reductions$umap@cell.embeddings 
    
    clist[[1]]@proj[[paste0('t-SNE', dim0)]] = my_tsne_proj
    clist[[1]]@proj[[paste0('UMAP', dim0)]] = my_umap_proj
  }
  
  
  saveRDS(eset, paste0(celloInputPath, '/eset.rds'))
  saveRDS(clist, paste0(celloInputPath, '/clist.rds'))
  if(!is.null(clusterOnDim) & !is.null(addPcaDims)) saveRDS(seurat.obj, seuratObjPath)
}



## update mtx of cello objects, without touch other stuff
## need provide feature names
updateRow4Cello <- function(celloInputPath, seuratObjPath, 
                            features_name, assay = 'ATAC'){
  eset = readRDS(paste0(celloInputPath, '/eset.rds'))
  mdata = phenoData(eset)
  
  seurat.obj = readRDS(seuratObjPath)
  DefaultAssay(seurat.obj) <- assay
  
  mtx = seurat.obj@assys[[assay]]@counts[features_name, ]
  norm_mtx = seurat.obj@assys[[assay]]@data[features_name, ]
  
  meta = seurat.obj@meta.data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx, 
                                       norm_exprs = norm_mtx),
              phenoData =  new("AnnotatedDataFrame", data = meta),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  
  phenoData(eset) <- mdata
  
  saveRDS(eset, paste0(celloInputPath, '/eset.rds'))
 }


## rbind new mtx to cello
rbindNewMtx4Cello <- function(celloInputPath, add_mtx, add_norm_mtx){
  eset = readRDS(paste0(celloInputPath, '/eset.rds'))
  mdata = phenoData(eset)
  mtx0 = eset@assayData$exprs
  norm_mtx0 = eset@assayData$norm_exprs
  
  add_mtx = add_mtx[, colnames(mtx0)]
  add_norm_mtx = add_norm_mtx[, colnames(mtx0)]
  
  mtx = rbind(mtx0, add_mtx)
  norm_mtx = rbind(norm_mtx0, add_norm_mtx)
  
  meta = mdata@data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx, 
                                       norm_exprs = norm_mtx),
              phenoData =  new("AnnotatedDataFrame", data = meta),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  
  phenoData(eset) <- mdata
  
  saveRDS(eset, paste0(celloInputPath, '/eset.rds'))
}


do_DA <- function(mtx_score, clusters, test = 'wilcox', fdr = 0.05, topn = 10){
  clusters$cluster = as.character(clusters$cluster)
  cls = unique(clusters$cluster)
  res = NULL
  features = rownames(mtx_score)
  for(cluster0 in cls){
    bc0 = clusters[cluster == cluster0]$barcode
    mtx1 = mtx_score[, colnames(mtx_score) %in% bc0]
    mtx2 = mtx_score[, !colnames(mtx_score) %in% bc0]
    mu1 = sapply(1:length(features), function(x) mean(mtx1[x, ]))
    mu2 = sapply(1:length(features), function(x) mean(mtx2[x, ]))
    
    
    pvs = sapply(1:length(features), function(x) wilcox.test(mtx1[x, ], mtx2[x, ], 
                                                             alternative = 'greater')$p.value )
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster' = cluster0,
                      'mean1' = mu1, 'mean2' = mu2,
                      'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}


#fg_genes: vector of forground genes
#bg_genes: background genes
#type: BP, CC, kegg
do_GO <- function(fg_genes, bg_genes = NULL, type = "BP", qCutoff = 0.05,
                  organism = c("mmu",  "hsa")) {
  if(organism =="mmu") {
    orgdb <- "org.Mm.eg.db"
    fromType = "SYMBOL"
    if(!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
    
  } else if(organism == "hsa") {
    orgdb <- "org.Hs.eg.db"
    if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
    fromType = "SYMBOL"
    
  }
  
  if(!is.null(bg_genes)) bg.df <- bitr(bg_genes, fromType = fromType,
                toType = c("SYMBOL", "ENTREZID"),
                OrgDb = orgdb)
  gene.df <- bitr(fg_genes, fromType = fromType,
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = orgdb)
  
  
  if(type == "kegg") {
    kegg_gene <- gene.df$ENTREZID
    if(!is.null(bg_genes)){
      kegg_bg <- bg.df$ENTREZID
      enrich_list <- enrichKEGG(
        gene          = kegg_gene,
        universe      = kegg_bg,
        organism      = organism,
        pAdjustMethod = "BH",
        qvalueCutoff  = qCutoff)
      
    }else{
      enrich_list <- enrichKEGG(
        gene          = kegg_gene,
        organism      = organism,
        pAdjustMethod = "BH",
        qvalueCutoff  = qCutoff)
      
    }
    
  }else {
    if(!is.null(bg_genes)){
      enrich_list <- enrichGO(gene        = gene.df$ENTREZID,
                              universe      = bg.df$ENTREZID,
                              OrgDb         = orgdb,
                              ont           = type,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = qCutoff,
                              readable      = TRUE)
    }else{
      enrich_list <- enrichGO(gene        = gene.df$ENTREZID,
                              OrgDb         = orgdb,
                              ont           = type,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = qCutoff,
                              readable      = TRUE)
    }
    
  }
  
  #return(enrich_list@result[enrich_list@result$qvalue <= qCutoff, ])
  return(enrich_list)
}

  
generate_gene_cisActivity <- function(gene_gtf, mtx, include_body = T,
                                      dist_to_tss = 2000){
  ## generating gene cis activity score
  ## input: gene gtf file, and atac matrix file
  
  ## 1. select gene up/down-stream regions (promoter with/without gene_body) ##
  
  gene_ann = fread(gene_gtf, sep = '\t')
  gene_ann = gene_ann[V3 == 'gene']
  gene_ann[, 'gene_name' := unlist(strsplit(V9, ';'))[3], by = V9]
  gene_ann[, 'gene_name' := gsub("\"", "", gene_name), by = gene_name]
  gene_ann[, 'gene_name' := unlist(strsplit(gene_name, ' '))[3], by = gene_name]
  names(gene_ann)[1] = 'chr'
  gene_ann = subset(gene_ann, select = c(chr, V4, V5, V7, gene_name))
  chrs = 1:22
  chrs = c(chrs, 'X', 'Y')
  gene_ann = gene_ann[chr %in% chrs]
  gene_ann = gene_ann[!duplicated(gene_name)]
  
  
  names(gene_ann)[2:4] = c('ss', 'ee', 'strand')
  if(!include_body){
    gene_ann[, 'start' := ifelse(strand == '+', ss - dist_to_tss, ee - dist_to_tss)]
    gene_ann[, 'end' := ifelse(strand == '+', ss + dist_to_tss, ee + dist_to_tss)]
  }else{
    gene_ann[, 'start' := ifelse(strand == '+', ss - dist_to_tss, ss)]
    gene_ann[, 'end' := ifelse(strand == '+', ee, ee + dist_to_tss)]
    
  }
  
  gene_ann[, 'chr' := paste0('chr', chr)]
  
  gene_ann = subset(gene_ann, select = c(chr, start, end, gene_name))
  
  
  ## 2. read mtx file ##
  
  rnames = rownames(mtx)
  chrs = sapply(rnames, function(x) unlist(strsplit(x, '-'))[1])
  starts = sapply(rnames, function(x) unlist(strsplit(x, '-'))[2])
  ends = sapply(rnames, function(x) unlist(strsplit(x, '-'))[3])
  
  peaks = data.table('chr' = chrs, 'start' = as.integer(starts), 
                     'end' = as.integer(ends))
  setkey(peaks, chr, start, end)
  peaks[, 'pname' := paste(chr, start, end, sep = '-')]
  over.ids = foverlaps(gene_ann, peaks, by.x = c('chr', 'start', 'end'),
                       by.y = c('chr', 'start', 'end'), which = T)
  over.ids[, 'gene_name' := gene_ann[xid, gene_name]]
  over.ids[, 'pname' := peaks[yid, pname]]
  over.ids = over.ids[complete.cases(over.ids)]
  
  
  smtx = sparseMatrix(i = over.ids$xid, j = over.ids$yid,
                      dimnames = list(gene_ann$gene_name[1:max(over.ids$xid)],
                                      peaks$pname[1:max(over.ids$yid)]))
  
  mtx = mtx[rownames(mtx) %in% colnames(smtx), ]
  smtx = smtx[, rownames(mtx)]
  
  activity.matrix = smtx %*% mtx
  rs = Matrix::rowSums(activity.matrix)
  activity.matrix = activity.matrix[rs > 10, ]
  
  return(activity.matrix)
  
}

## calculat gene activity score given interactions and atac_mtx
## conn = data.table(Peak1, Peak2, coaccess)
gene_activity_gConn <- function(conn, gene_ann, atac_mtx,
                                distal_dist = 2e5){
  ## rm connections that are too far from each other
  ## and construct p2p matrix for idication of interacton
  conns1 <- tidyr::separate(conn, col = Peak1, 
                            into = c('chr1', 'start1', 'end1'))
  conns2 <- tidyr::separate(conn, col = Peak2, 
                            into = c('chr2', 'start2', 'end2'))
  conns1$coaccess = NULL
  conn <- cbind(conns1, conns2)
  class(conn$start1) = class(conn$end1) = 
    class(conn$start2) = class(conn$end2) = 'integer'
  conn = conn[chr1 == chr2]
  conn = conn[abs(start2/2 + end1/2 - start2/2 - end2/2) <= distal_dist]
  
  upeaks = sort(unique(union(unique(conn$Peak2), unique(conn$Peak1))))
  id1 = sapply(conn$Peak1, function(x) which(upeaks == x))
  id2 = sapply(conn$Peak2, function(x) which(upeaks == x))
  
  p2p.mtx <- sparseMatrix(i = id1, j = id2,
                          dimnames = list(upeaks, upeaks))
  p2p.mtx = p2p.mtx * 1
  
  ## construct gene peak overlapping idication matrix
  peaks = rownames(atac_mtx)
  peaks = sapply(peaks, function(x) unlist(strsplit(x, ','))[1])
  rownames(atac_mtx) <- peaks
  
  peaks = peaks[peaks %in% upeaks]
  peaks.dt = tidyr::separate(data.table(x = peaks), col = x,
                             into = c('chr', 'start', 'end'))
  peaks.dt$peak_name = peaks
  class(peaks.dt$start) = 'integer'
  class(peaks.dt$end) = 'integer'
  setkey(peaks.dt, chr, start, end)
  
  
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  gene_ann[, 'start' := tss - 1000]
  gene_ann[, 'end' := tss + 1000]
  over.ids = foverlaps(gene_ann, peaks.dt, by.x = c('chr', 'start', 'end'),
                       by.y = c('chr', 'start', 'end'), which = T)
  over.ids[, 'gene_name' := gene_ann[xid, gene_name]]
  over.ids[, 'peak_name' := peaks.dt[yid, peak_name]]
  over.ids = over.ids[complete.cases(over.ids)]
  
  smtx = sparseMatrix(i = over.ids$xid, j = over.ids$yid,
                      dimnames = list(gene_ann$gene_name[1:max(over.ids$xid)],
                                      peaks.dt$peak_name[1:max(over.ids$yid)]))
  
  ## sum the same gene
  rn = rownames(smtx)
  dup.rn = unique(rn[duplicated(rn)])
  uniq.rn = setdiff(unique(rn), dup.rn)
  smtx.u <- smtx[uniq.rn, ]
  
  smtx.d <- sapply(dup.rn, function(x){
    tmpID = which(rn == x)
    return(Matrix::colSums(smtx[tmpID, ]))
  })
  
  g2p.mtx = rbind(smtx.u, t(smtx.d) > 0)
  
  ## get the score
  p2p.mtx = p2p.mtx[colnames(g2p.mtx), ]
  p2p.mtx = p2p.mtx[, colnames(g2p.mtx)]
  atac_mtx = atac_mtx[colnames(g2p.mtx), ]
  gascore.mtx = g2p.mtx %*% p2p.mtx %*% atac_mtx
  return(gascore.mtx)
}


## informative genes selection using gini index from qin
## if gini_cut_qt is an integer, it will output top 1:gini_cut_qt genes as include genes
ifg_select <- function(data, cluster, cluster_min_cell_num = 100, 
                       min_cluster_expr_fraction = .1, gini_cut_qt = .75, 
                       do_go = F, filePath = NULL, 
                       fileName = "ifg_select_",
                       orgdb = "org.Hs.eg.db", gene_id_type = "SYMBOL",
                       gene_background = NULL) {
  
  if(ncol(data) != length(cluster)) stop("Cell number do not match cluster length.")
  
  use_clus <- names(which(table(cluster) >= cluster_min_cell_num))
  
  # For each cluster, compute gene expressed fraction
  
  expr_clus_frac <- sapply(use_clus, function(x) {
    
    cur_data <- data[, cluster == x]
    
    Matrix::rowMeans(cur_data > 0)
    
  })
  
  # Compute gini coefficient 
  
  # Require a gene to be expressed in at least one cluster with at least .1 expressed fraction to be considered for downstream uniform gene selection
  
  use_g <- rownames(expr_clus_frac)[rowSums(expr_clus_frac >= min_cluster_expr_fraction) > 0] # 11813
  
  message(paste0("Selecting informative features from ", length(use_g), " robustly detected features."))
  
  expr_clus_frac <- expr_clus_frac[rownames(expr_clus_frac) %in% use_g,]
  
  gene_clus_gini <- apply(expr_clus_frac, 1, gini)
  
  
  
  #pdf(paste0(filePath, fileName, "gene_clus_gini_hist.pdf"))
  
  #hist(gene_clus_gini, breaks = 100)
  
  #dev.off()
  
  
  
  gene_clus_gini = sort(gene_clus_gini, decreasing = T)
  if(gini_cut_qt > 1){
    topn = min(length(use_g), gini_cut_qt)
    exclude_g <- names(gene_clus_gini[-(1:topn)])
    
    include_g <- names(gene_clus_gini[(1:topn)])
  }else{
    gini_cut <- quantile(gene_clus_gini, gini_cut_qt) 
    
    message(paste0("Cut at gini quantile ", gini_cut_qt, " with value ", gini_cut))
    
    exclude_g <- names(gene_clus_gini)[gene_clus_gini < gini_cut]
    
    include_g <- names(gene_clus_gini)[gene_clus_gini >= gini_cut]
    
  }
  
  #write.csv(data.frame(less_specific_feature = exclude_g), paste0(filePath, fileName, "less_specific_feature_list.csv"))
  
  #write.csv(data.frame(specific_feature = include_g), paste0(filePath, fileName, "specific_feature_list.csv"))
  
  
  
  message(paste0("Found ", length(exclude_g), " less specific features."))
  
  message(paste0("Returning ", length(include_g), " specific features."))
  
  
  
  if(do_go) {
    
    if(!length(gene_background)) {
      
      gene_background <- use_g
      
    }
    
    exclude_g.df <- bitr(exclude_g, fromType = gene_id_type,
                         
                         toType = c("ENTREZID"),
                         
                         OrgDb = orgdb)
    
    bg.df <- bitr(gene_background, fromType = gene_id_type,
                  
                  toType = c("ENTREZID"),
                  
                  OrgDb = orgdb)
    
    exclude_g_go<- enrichGO(gene        = exclude_g.df$ENTREZID,
                            
                            universe      = bg.df$ENTREZID,
                            
                            OrgDb         = orgdb,
                            
                            ont           = "BP",
                            
                            pAdjustMethod = "BH",
                            
                            pvalueCutoff  = 0.01,
                            
                            qvalueCutoff  = 0.05,
                            
                            readable      = TRUE)
    
    write.csv(exclude_g_go, paste0(filePath, fileName, "exlude_feature_go.csv"))
    
    
    
    
    
    include_g.df <- bitr(include_g, fromType = gene_id_type,
                         
                         toType = c("ENTREZID"),
                         
                         OrgDb = orgdb)
    
    bg.df <- bitr(gene_background, fromType = gene_id_type,
                  
                  toType = c("ENTREZID"),
                  
                  OrgDb = orgdb)
    
    include_g_go <- enrichGO(gene        = include_g.df$ENTREZID,
                            
                            universe      = bg.df$ENTREZID,
                            
                            OrgDb         = orgdb,
                            
                            ont           = "BP",
                            
                            pAdjustMethod = "BH",
                            
                            pvalueCutoff  = 0.01,
                            
                            qvalueCutoff  = 0.05,
                            
                            readable      = TRUE)
    
    write.csv(include_g_go, paste0(filePath, fileName, "include_feature_go.csv"))
    
  }
  
  
  
  return(list('include_g' = include_g, 'exclude_g' = exclude_g))
  
}


## re-select genes/features by F-stat
fstat_select <- function(data, cluster, cluster_min_cell_num = 100, 
                       min_cluster_expr_fraction = .1, cut_qt = 0.75, 
                       do_go = F, filePath = NULL, 
                       fileName = "fstat_select_",
                       orgdb = "org.Hs.eg.db", gene_id_type = "SYMBOL",
                       gene_background = NULL,
                       max_cell_per_cl = 200) {
  
  if(ncol(data) != length(cluster)) stop("Cell number do not match cluster length.")
  nn = names(cluster)
  cluster = as.character(cluster)
  names(cluster) = nn
  
  use_clus <- names(which(table(cluster) >= cluster_min_cell_num))
  
  # For each cluster, compute gene expressed fraction
  
  expr_clus_frac <- sapply(use_clus, function(x) {
    
    cur_data <- data[, cluster == x]
    
    Matrix::rowMeans(cur_data > 0)
    
  })
  
  # Compute F-stats
  
  # Require a gene to be expressed in at least one cluster with at least .1 expressed fraction to be considered for downstream uniform gene selection
  
  use_g <- rownames(expr_clus_frac)[rowSums(expr_clus_frac >= min_cluster_expr_fraction) > 0] 
  
  message(paste0("Selecting informative features from ", length(use_g), " robustly detected features."))
  
  expr_clus_frac <- expr_clus_frac[rownames(expr_clus_frac) %in% use_g,]
  
  data = data[rownames(data) %in% use_g, ]
  data = log1p(data) 
  
  cluster = cluster[cluster %in% use_clus]
  data = data[, names(cluster)]
  
  ##downsample each cluster to have at most 200 cells
  cluster0 = NULL
  set.seed(2019)
  for(cl in use_clus){
    tmp = cluster[cluster == cl]
    if(length(tmp) > max_cell_per_cl){
      cluster0 = c(cluster0, tmp[sample(1:length(tmp), max_cell_per_cl)])
    }else{
      cluster0 = c(cluster0, tmp)
    }
  }
  
  data = data[, names(cluster0)]
  data = (data > 0) * 1
  gene_clus_fstat <- apply(data, 1, function(x) anova(lm(x ~ cluster0))$`F value`[1])
  #gene_clus_fstat = rep(0, nrow(data))
  #for(i in 1:nrow(data)){
  #  gene_clus_fstat[i] = anova(lm(data[i, ] ~ cluster))$`F value`[1]
  #}
  
  #gene_clus_fstat = pmin(gene_clus_fstat * length(gene_clus_fstat), 1)
  
  pdf(paste0(filePath, fileName, "gene_clus_fstat_pv_hist.pdf"))
  
  hist(gene_clus_fstat, breaks = 100)
  
  dev.off()
  
  cut_thr = quantile(gene_clus_fstat, cut_qt)
  exclude_g <- names(gene_clus_fstat)[gene_clus_fstat < cut_thr]
  
  include_g <- names(gene_clus_fstat)[gene_clus_fstat >= cut_thr]
  
  #write.csv(data.frame(less_specific_feature = exclude_g), paste0(filePath, fileName, "less_specific_feature_list.csv"))
  
  #write.csv(data.frame(specific_feature = include_g), paste0(filePath, fileName, "specific_feature_list.csv"))
  
  
  
  message(paste0("Found ", length(exclude_g), " less specific features."))
  
  message(paste0("Returning ", length(include_g), " specific features."))
  
  
  
  if(do_go) {
    
    if(!length(gene_background)) {
      
      gene_background <- use_g
      
    }
    
    exclude_g.df <- bitr(exclude_g, fromType = gene_id_type,
                         
                         toType = c("ENTREZID"),
                         
                         OrgDb = orgdb)
    
    bg.df <- bitr(gene_background, fromType = gene_id_type,
                  
                  toType = c("ENTREZID"),
                  
                  OrgDb = orgdb)
    
    exclude_g_go<- enrichGO(gene        = exclude_g.df$ENTREZID,
                            
                            universe      = bg.df$ENTREZID,
                            
                            OrgDb         = orgdb,
                            
                            ont           = "BP",
                            
                            pAdjustMethod = "BH",
                            
                            pvalueCutoff  = 0.01,
                            
                            qvalueCutoff  = 0.05,
                            
                            readable      = TRUE)
    
    write.csv(exclude_g_go, paste0(filePath, fileName, "exlude_feature_go.csv"))
    
    
    
    
    
    include_g.df <- bitr(include_g, fromType = gene_id_type,
                         
                         toType = c("ENTREZID"),
                         
                         OrgDb = orgdb)
    
    bg.df <- bitr(gene_background, fromType = gene_id_type,
                  
                  toType = c("ENTREZID"),
                  
                  OrgDb = orgdb)
    
    include_g_go <- enrichGO(gene        = include_g.df$ENTREZID,
                             
                             universe      = bg.df$ENTREZID,
                             
                             OrgDb         = orgdb,
                             
                             ont           = "BP",
                             
                             pAdjustMethod = "BH",
                             
                             pvalueCutoff  = 0.01,
                             
                             qvalueCutoff  = 0.05,
                             
                             readable      = TRUE)
    
    write.csv(include_g_go, paste0(filePath, fileName, "include_feature_go.csv"))
    
  }
  
  
  
  return(list('include_g' = include_g, 'exclude_g' = exclude_g))
  
}




#integrate a list of seurat objects
#supporse each element of the list was normalized and/or scaled
integrateSeuratList <- function(seurat.list, npc = 50, reg.var = NULL,
                               reduction = 'pca', anchor.features = 6000){
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                           anchor.features = anchor.features)
  rm(seurat.list)
  
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:npc)
  rm(seurat.anchors)
  
  DefaultAssay(seurat.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE, vars.to.regress = reg.var)
  seurat.integrated <- RunPCA(seurat.integrated, npcs = npc, verbose = FALSE)
  #seurat.integrated <- RunTSNE(seurat.integrated, reduction = reduction, dims = 1:npc,
  #                             check_duplicates = FALSE)
  #seurat.integrated <- RunUMAP(seurat.integrated, reduction = reduction, dims = 1:npc)
 
  return(seurat.integrated)
  
}

#integrate a list of seurat objects
integrateSeuratList_withReference <- function(seurat.list, npc = 50, reg.var = NULL,
                                reduction = 'pca', anchor.features = 6000,
                                ref.name = 'sample1'){
  ref.dataset <- which(names(seurat.list) == ref.name)
  
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                           anchor.features = anchor.features,
                                           reference = ref.dataset)
  rm(seurat.list)
  
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:npc)
  rm(seurat.anchors)
  
  DefaultAssay(seurat.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE, vars.to.regress = reg.var)
  seurat.integrated <- RunPCA(seurat.integrated, npcs = npc, verbose = FALSE)
  #seurat.integrated <- RunTSNE(seurat.integrated, reduction = reduction, dims = 1:npc,
  #                             check_duplicates = FALSE, 
  #                             tsne.method = "FIt-SNE", nthreads = 2, max_iter = 2000)
  #seurat.integrated <- RunUMAP(seurat.integrated, reduction = reduction, dims = 1:npc)
  
  return(seurat.integrated)
  
}



#integrate a list of seurat objects by SCT normalization
integrateSeuratList_SCT <- function(seurat.list, npc = 50, reg.var = NULL,
                                reduction = 'pca', anchor.features = 3000,
                                resolution = 0.6){
  anchor.features = SelectIntegrationFeatures(seurat.list, 
                                            nfeatures = anchor.features)
  seurat.list <- PrepSCTIntegration(object.list = seurat.list, 
                                           verbose = F,
                                           anchor.features = anchor.features)
  seurat.anchors = FindIntegrationAnchors(object.list = seurat.list, 
                                          normalization.method = "SCT", 
                                          anchor.features = anchor.features, 
                                          verbose = FALSE)
  rm(seurat.list)
  
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:npc,
                                     normalization.method = "SCT", verbose = F)
  rm(seurat.anchors)
  
  DefaultAssay(seurat.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  seurat.integrated <- RunPCA(seurat.integrated, npcs = npc, verbose = FALSE)
  seurat.integrated <- RunTSNE(seurat.integrated, reduction = reduction, dims = 1:npc)
  seurat.integrated <- RunUMAP(seurat.integrated, reduction = reduction, dims = 1:npc)
  seurat.integrated <- FindNeighbors(seurat.integrated, reduction = reduction, 
                                     dims = 1:npc)
  seurat.integrated <- FindClusters(seurat.integrated, resolution = resolution)
  
  return(seurat.integrated)
  
}


## define ctype signature score by marker genes
ctype_assign_chisq <- function(mtx, ctype_markers, 
                               upper_thr = 0.1, lower_thr = 0.5,
                               negative_marker_w = 0.3,
                               doSmooth = F, cell_embedding = NULL,
                               k = 20){
  ctypes = unique(ctype_markers$ctype)
  
  mtx = as.matrix(mtx[rownames(mtx) %in% unique(ctype_markers$gene), ])
  score.mtx = chisq.test(mtx)$residuals
  score.mtx[is.na(score.mtx)] = 0
  
  score.mtx = as.matrix(score.mtx[rownames(score.mtx) %in% unique(ctype_markers$gene), ])
  score.mtx.neg = score.mtx.pos = score.mtx
  for(i in 1:nrow(score.mtx)){
    tmp = score.mtx[i,]
    uthr = quantile(tmp[tmp > 0], upper_thr)
    tmp_pos = ifelse(tmp >= uthr, 1, 0)
    dthr = quantile(tmp[tmp < 0], lower_thr)
    tmp_neg = ifelse(tmp <= dthr, 1, 0)
    
    score.mtx.pos[i, ] = tmp_pos
    score.mtx.neg[i, ] = tmp_neg
    
    
  }
  
  # calculate support freq for each marker gene given a ctype
  supp_freq = list()
  score_cell = list()
  for(ctype0 in ctypes){
    genes0_freq = ctype_markers[ctype == ctype0, ]
    genes0_freq = genes0_freq[order(-N), ]
    if(length(which(rownames(score.mtx) %in% genes0_freq$gene)) == 0) next
    genes0_freq = genes0_freq[gene %in% rownames(score.mtx), ]
    
    # different weight for neg and postive markers
    genes0_freq[, freq := ifelse(direction == -1, 
                                 freq * negative_marker_w, freq *(1-negative_marker_w))]
    genes0_freq$freq = genes0_freq$freq/sum(genes0_freq$freq)
    genes0_freq.pos = genes0_freq[direction == 1]
    genes0_freq.neg = genes0_freq[direction == -1]
    
    
    score.mtx0.pos = score.mtx.pos[rownames(score.mtx.pos) %in% genes0_freq.pos$gene, ]
    
    if(nrow(genes0_freq.pos) == 1){
      score_cell[[ctype0]] = score.mtx0.pos * genes0_freq.pos$freq 
    }else{
      score_cell[[ctype0]] = t(genes0_freq.pos$freq) %*% score.mtx0.pos
    }
    
    if(any(genes0_freq$direction == -1)){
      score.mtx0.neg = score.mtx.neg[rownames(score.mtx.pos) %in% genes0_freq.neg$gene, ]
      
      if(nrow(genes0_freq.neg) == 1){
        score_cell[[ctype0]] = score_cell[[ctype0]] + score.mtx0.neg * genes0_freq.neg$freq 
      }else{
        score_cell[[ctype0]] = score_cell[[ctype0]] + t(genes0_freq.neg$freq) %*% score.mtx0.neg
      }
    }
    
    
    
    
  }
  
  res = do.call('rbind', score_cell)
  rownames(res) = names(score_cell)
  
  ##smooth by cell_embedding
  if(doSmooth){
    if(is.null(cell_embedding)) stop('cell embedding mtx needed for smoothing!')
    knn.ids <- FNN::knn.index(cell_embedding, k = k)
    mtx.ids <- sparseMatrix(i = rep(c(1:ncol(res)), each = k) , 
                            j = as.vector(t(knn.ids)),
                            dimnames = list(colnames(res), 
                                            colnames(res)))
    mtx.ids = mtx.ids + diag(ncol(mtx.ids))
    res = res %*% mtx.ids/(k+1)
  }
  res[res > 1] = 1
  return(res)
}


## define ctype signature score by marker genes --stringent version
ctype_assign_chisq_update <- function(mtx, ctype_markers, 
                               upper_thr = 0.1, lower_thr = 0.1,
                               negative_marker_w = 0.5){
  ctypes = unique(ctype_markers$ctype)
  
  mtx = as.matrix(mtx[rownames(mtx) %in% unique(ctype_markers$gene), ])
  score.mtx = chisq.test(mtx)$residuals
  score.mtx[is.na(score.mtx)] = 0
  
  score.mtx = as.matrix(score.mtx[rownames(score.mtx) %in% unique(ctype_markers$gene), ])
  score.mtx.neg = score.mtx.pos = score.mtx
  for(i in 1:nrow(score.mtx)){
    tmp = score.mtx[i,]
    uthr = quantile(tmp[tmp > 0], upper_thr)
    tmp_pos = ifelse(tmp >= uthr, 1, 0)
    dthr = quantile(tmp[tmp > 0], lower_thr)
    tmp_neg = ifelse(tmp <= dthr, 1, 0)
    
    score.mtx.pos[i, ] = tmp_pos
    score.mtx.neg[i, ] = tmp_neg
    
    
  }
  
  # calculate support freq for each marker gene given a ctype
  supp_freq = list()
  score_cell = list()
  for(ctype0 in ctypes){
    genes0_freq = ctype_markers[ctype == ctype0, ]
    genes0_freq = genes0_freq[order(-N), ]
    if(length(which(rownames(score.mtx) %in% genes0_freq$gene)) == 0) next
    genes0_freq = genes0_freq[gene %in% rownames(score.mtx), ]
    
    # different weight for neg and postive markers
    genes0_freq[, freq := ifelse(direction == -1, 
                                 freq * negative_marker_w, freq *(1-negative_marker_w))]
    genes0_freq$freq = genes0_freq$freq/sum(genes0_freq$freq)
    genes0_freq.pos = genes0_freq[direction == 1]
    genes0_freq.neg = genes0_freq[direction == -1]
    
    
    score.mtx0.pos = score.mtx.pos[rownames(score.mtx.pos) %in% genes0_freq.pos$gene, ]
    
    if(nrow(genes0_freq.pos) == 1){
      score_cell[[ctype0]] = score.mtx0.pos * genes0_freq.pos$freq 
    }else{
      score_cell[[ctype0]] = t(genes0_freq.pos$freq) %*% score.mtx0.pos
    }
    
    if(any(genes0_freq$direction == -1)){
      score.mtx0.neg = score.mtx.neg[rownames(score.mtx.pos) %in% genes0_freq.neg$gene, ]
      
      if(nrow(genes0_freq.neg) == 1){
        score_cell[[ctype0]] = score_cell[[ctype0]] + score.mtx0.neg * genes0_freq.neg$freq 
      }else{
        score_cell[[ctype0]] = score_cell[[ctype0]] + t(genes0_freq.neg$freq) %*% score.mtx0.neg
      }
    }
    
    
    
    
  }
  
  res = do.call('rbind', score_cell)
  rownames(res) = names(score_cell)
  return(res)
}






## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, myColors,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) + scale_fill_manual(values = myColors)
  return(p)
}


## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, myColors,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, myColors = myColors,
                                                              pt.size = pt.size, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



signatureScore_zscore <- function(seurat.obj, pfeatures, nfeatures = NULL,
                                  vars.to.regress = NULL,
                                  score.name = 'quicent'){
  
  pfeatures = pfeatures[pfeatures %in% rownames(seurat.obj)]
  if(!is.null(nfeatures)) nfeatures = nfeatures[nfeatures %in% rownames(seurat.obj)]
  
  seurat.obj <- ScaleData(seurat.obj, features = c(pfeatures, nfeatures),
                            do.scale = T, do.center = T, 
                          vars.to.regress = vars.to.regress)
  
  mtx = GetAssayData(seurat.obj, slot = 'scale.data')
  
  pscore = mtx[pfeatures, ,drop = F]
  
  if(length(nfeatures) > 0) nscore = mtx[nfeatures, ,drop = F]
  
  pscore = Matrix::colSums(pscore)
  
  if(length(nfeatures) > 0){
    nscore = Matrix::colSums(nscore)
    scores = (pscore - nscore)/length(c(pfeatures, nfeatures))
  }else{
    scores = pscore/length(pfeatures)
  }
  seurat.obj <- AddMetaData(seurat.obj, metadata = scores, 
                            col.name = score.name)
  return(seurat.obj)
}

## use robust zscore, (x-median)/mad
signatureScore_zscore_robust <- function(seurat.obj, pfeatures, nfeatures = NULL,
                                  vars.to.regress = NULL,
                                  score.name = 'quicent'){
  
  pfeatures = pfeatures[pfeatures %in% rownames(seurat.obj)]
  if(!is.null(nfeatures)) nfeatures = nfeatures[nfeatures %in% rownames(seurat.obj)]
  
  ## regress out
   seurat.obj <- ScaleData(seurat.obj, features = c(pfeatures, nfeatures),
                          do.scale = F, do.center = F, 
                          vars.to.regress = vars.to.regress)
  
  mtx = GetAssayData(seurat.obj, slot = 'scale.data')
  
  ## calculated robust zscore
  rzscore <- function(x){
    sd.robust = sd(x)/2 + mad(x)/2
    x.robust = (x - median(x))/sd.robust
    if(sd.robust <= 0) x.robust = x - median(x)
    return(x.robust)
  }
  
  pscore = mtx[pfeatures, ,drop = F]
  pscore = t(apply(pscore, 1, rzscore))
  pscore = Matrix::colSums(pscore)
  
  if(length(nfeatures) > 0) {
    nscore = mtx[nfeatures, ,drop = F]
    nscore = t(apply(nscore, 1, rzscore))
    nscore = Matrix::colSums(nscore)
    scores = (pscore - nscore)/length(c(pfeatures, nfeatures))
  }else{
    scores = pscore/length(pfeatures)
  }
  seurat.obj <- AddMetaData(seurat.obj, metadata = scores, 
                            col.name = score.name)
  return(seurat.obj)
}

## remove features active in less than min_frac_per_cluser in first class
## do downsample to max_cellPer_cluster
## support wilcox and t test only
## test-model: one-rest; control
runDiffMotifEnrich <- function(mtx_score, clusters, test = 'wilcox',
                        fdr = 0.01, topn = 20,
                        min_frac_per_cluster = 0.2,
                        max_cell_per_clust = 300,
                        test.mode = 'one-rest', control = NULL){
  set.seed(2020)
  clusters$cluster = as.character(clusters$cluster)
  cls = unique(clusters$cluster)
  res = NULL
  features = rownames(mtx_score)
  mtx_score = mtx_score[, clusters$barcode]
  for(cluster0 in cls){
    bc0 = clusters[cluster == cluster0]$barcode
    mtx1 = mtx_score[, colnames(mtx_score) %in% bc0]
    if(test.mode == 'one-rest') {
      mtx2 = mtx_score[, !colnames(mtx_score) %in% bc0]
    }
    if(test.mode == 'control'){
      if(is.null(control)) stop('Should specific control cluster!')
      bc2 = clusters[cluster %in% control]$barcode
      mtx2 = mtx_score[, colnames(mtx_score) %in% bc2]
      if(cluster0 %in% control) next
    }
    mu1 = sapply(1:length(features), function(x) mean(mtx1[x, ]))
    mu2 = sapply(1:length(features), function(x) mean(mtx2[x, ]))
    
    pvs = rep(0.5, length(features))
    
    for(x in 1:length(features)){
      a1 = mtx1[x, ]
      a2 = mtx2[x, ]
      if(mu1[x] <= 0) next
      #if(mu1[x] <= mu2[x]) next
      if(length(which(!is.na(a1))) < 2 || length(which(!is.na(a2))) < 2) next
      if(mean(a1 > 0) < min_frac_per_cluster) next
      if(!is.null(max_cell_per_clust)){
        if(length(a1) > max_cell_per_clust) a1 = sample(a1, max_cell_per_clust)
        if(length(a2) > max_cell_per_clust) a2 = sample(a2, max_cell_per_clust)
      }
      
      if(test == 'wilcox') pvs[x] = wilcox.test(a1, a2, alternative = 'greater')$p.value
      if(test == 't') pvs[x] = t.test(a1, a2, alternative = 'greater')$p.value
      
    }
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster1' = cluster0,
                      'mean1' = mu1, 'mean0' = mu2, 
                      'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    res0 = res0[mean1 > 0]
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}


## mtx_objs: a list of matrix 
run_integration <- function(mtx_list, integrate_by = 'VFACS',
                            top_variable_features = 5000, 
                            norm_by = 'tf-idf', nREDUCTION = 30,
                            minFrac_in_cell = 0.01, min_depth = 1000,
                            max_depth = 50000, reg.var = 'nCount_ATAC',
                            anchor.features = 2000, ref.sample4seurat = 1,
                            resolution = 0.6, verbose = F){
  nsample = length(mtx_list)
  sampleNames = names(mtx_list)
  if(is.null(sampleNames)) sampleNames = paste0('sample', 1:nsample)
  names(mtx_list) = sampleNames
  seu.all <- list()
  for(sample0 in sampleNames){
    # filter each mtx
    mtx_list[[sample0]] <- filterMat(mtx_list[[sample0]], minFrac_in_cell = minFrac_in_cell,
                                     min_depth = min_depth, max_depth = max_depth)
    if(integrate_by %in% c('seurat')){
      # create a seurat obj for each sample
      seurat.obj = runSeurat_Atac(mtx_list[[sample0]], npc = nREDUCTION, norm_by = norm_by, 
                                  top_variable_features = top_variable_features, 
                                  reg.var = reg.var)
      
      seurat.obj$sample = sample0
      
      seu.all[[sample0]] = seurat.obj
    }
  }
  
  ## pool/integrate data into a seurat object
  if(integrate_by == 'seurat'){
    refs = seu.all[ref.sample4seurat]
    seurat.obj <- FindIntegrationAnchors(object.list = seu.all[-ref.sample4seurat],
                                         reference = refs,
                                         anchor.features = anchor.features)
    rm(seu.all)
    seurat.obj <- IntegrateData(anchorset = seurat.obj, dims = 1:nREDUCTION)
    DefaultAssay(seurat.obj) <- "integrated"
    seurat.obj <- ScaleData(seurat.obj, verbose = FALSE,
                            features = VariableFeatures(seurat.obj))
    seurat.obj <- RunPCA(seurat.obj, npcs = nREDUCTION, verbose = verbose)
  }else{
    # pool the matrxi first with union features
    nf = sapply(mtx_list, nrow)
    nc = sapply(mtx_list, ncol)
    
    umtx <- cBind_union_features(mtx_list)
    
    rm(mtx_list)
    seurat.obj = runSeurat_Atac(umtx, npc = nREDUCTION, norm_by = norm_by, 
                                top_variable_features = top_variable_features, 
                                reg.var = reg.var)
    seurat.obj$sample = rep(sampleNames, nc)
    
  }
  
  if(integrate_by == 'pool') seurat.obj <- regress_on_pca(seurat.obj, 'sample')
  
  if(integrate_by == 'VFACS'){
    ## cluster and then reselect features
    ## variable features across clusters
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:nREDUCTION, reduction = 'pca', 
                                verbose = verbose)
    seurat.obj <- FindClusters(seurat.obj, resl = resolution, verbose = verbose)
    clusters = as.character(seurat.obj$seurat_clusters)
    mtx = seurat.obj@assays$ATAC@counts
    mtx_by_cls <- sapply(unique(clusters), function(x) {
      
      cl_data <- mtx[, clusters == x]
      
      Matrix::rowMeans(cl_data > 0)
      
    })
    mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
    sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
    names(sds) = rownames(mtx_by_cls.norm)
    sele.features = names(which(sds >= sort(sds, decreasing = T)[top_variable_features]))
    mtx0 = mtx[sele.features, ]
    mtx0.norm = TF_IDF(mtx0)
    seurat.obj@assays$ATAC@data[sele.features, ] <- mtx0.norm
    VariableFeatures(seurat.obj) <- sele.features
    seurat.obj <- RunPCA(seurat.obj, dims = 1:nReduction, verbose = verbose)
    seurat.obj <- regress_on_pca(seurat.obj, reg.var = reg.var)
    seurat.obj <- FindNeighbors(seurat.obj, verbose = verbose, 
                                dims = 1:nREDUCTION, reduction = 'pca')
    seurat.obj <- FindClusters(seurat.obj, verbose = verbose, resl = resolution)
    
  }
  
  seurat.obj <- RunUMAP(seurat.obj, reduction = "pca", 
                        verbose = verbose, dims = 1:nREDUCTION)
  
  
  if(integrate_by == 'harmony'){
    
    seurat.obj <- harmony::RunHarmony(seurat.obj, c("sample"), assay.use = 'ATAC')
    
    seurat.obj <- seurat.obj %>% 
      RunUMAP(reduction = "harmony", dims = 1:nREDUCTION, verbose = verbose) 
  }
  
  ## clustering on the integrated data
  ## seurat implemented louvain algorithm
  redm = ifelse(integrate_by == 'harmony', 'harmony', 'pca')
  seurat.obj = FindNeighbors(seurat.obj, reduction = redm, 
                             dims = 1:nREDUCTION, verbose = verbose)
  
  seurat.obj = FindClusters(seurat.obj, resolution = resolution, verbose = verbose)
  seurat.obj$active_clusters = seurat.obj$seurat_clusters
  
  
  return(seurat.obj)
}

## calculate person residual as in paper 
## 'Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data'
## By Dmitry Kobak, Genome Biology 2021
## given a gene-by-cell count matrix
pearson_residual <- function(mat, theta = 100){
  N = sum(mat)
  rs = rowSums(mat)
  cs = colSums(mat)
  mu = rs %o% cs/N
  z = (mat - mu)/sqrt(mu + mu^2/theta)
  rnull = names(which(rs == 0))
  z[rnull, ] <- 0
  z = Matrix(z, sparse = T)
  return(z)
}
