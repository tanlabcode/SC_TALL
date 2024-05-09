source('scDataAnalysis_Utilities_simp.R')
library(igraph)
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

library(Rcpp)
sourceCpp(code = '
          #include <Rcpp.h>
          using namespace Rcpp;
          // get overlap between two data frame, output a vector, in which 
          // the ith component to be 1 if the ith row of data1 has overlap in data2, zero otherwise
          // [[Rcpp::export]]
          IntegerVector getOverlaps_C1(DataFrame dat1, DataFrame dat2) {
          CharacterVector chr1 = dat1["chr"];
          CharacterVector chr2 = dat2["chr"];
          NumericVector start1 = dat1["start"];
          NumericVector end1 = dat1["end"];
          NumericVector start2 = dat2["start"];
          NumericVector end2 = dat2["end"];
          
          int n1 = chr1.size(), n2 = chr2.size();
          NumericVector midP1(n1), len1(n1), len2(n2), midP2(n2);
          IntegerVector over1(n1);
          
          len1 = (end1 - start1)/2;
          midP1 = (end1 + start1)/2;
          
          len2 = (end2 - start2)/2;
          midP2 = (end2 + start2)/2;
          
          for(int i=0; i<n1; i++){
          over1[i] = 0;
          for(int j=0; j<n2; j++){
          if((chr2[j] == chr1[i]) && (fabs(midP1[i] - midP2[j]) <= max(NumericVector::create(len1[i], len2[j])))){
          over1[i] = 1;
          break;
          }
          }
          }
          
          return(over1);
          }
')


## virids style, downsample, match_cell, colomn color by sample/cluster
plot_enrich_tf1 <- function(sele.tfs, zscore.all, bc_clusters,
                            cluster.levels = NULL, 
                            cluster.col = NULL,
                            up.qt = 0.95, low.qt = 0.05,
                            ndownsample = 2000, match_cell = F,
                            reScale = F, cluster.rows = F,
                            color_style = 'virid'){
  sele.zscores = zscore.all[sele.tfs, ]
  
  if(reScale) sele.zscores = t(scale(t(sele.zscores), center = T, scale = T))
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample and match_cell
  ncell.cl = min(table(bc_clusters$cluster))
  set.seed(2020)
  bc_clusters.down = bc_clusters
  if(match_cell){
    bc_clusters.down = NULL
    for(cl0 in unique(bc_clusters$cluster)){
      tmp = bc_clusters[cluster == cl0]
      if(nrow(tmp) > ncell.cl) tmp = tmp[sort(sample((1:nrow(tmp)), ncell.cl)), ]
      bc_clusters.down = rbind(bc_clusters.down, tmp)
    }
  }
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = bc_clusters$cluster,
                          'barcode' = bc_clusters$barcode,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  nc = length(unique(bc_clusters$cluster))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(is.null(cluster.col)){
    if(nc >= 3) color_cluster = getPalette(nc)
    if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
    names(color_cluster) = sort(unique(bc_clusters$cluster))
  }else{
    color_cluster = cluster.col
  }
  
  
  ## order 
  ann_column = ann_column[order(factor(ann_column$cluster,
                                       levels = sort(unique(ann_column$cluster)))), ]
  if(is.null(cluster.levels)) {
    
    ann_column = ann_column[order(factor(ann_column$cluster,
                                         levels = sort(unique(ann_column$cluster)))), ]
    
  }else{
    ann_column = ann_column[order(factor(ann_column$cluster,
                                         levels = cluster.levels)), ]
    
    color_cluster = color_cluster[cluster.levels]
  }
  
  
  ann_colors = list('cluster' = color_cluster)
  
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode <- NULL
  
  color_fun = viridis(100)
  if(color_style == 'purple-yellow') color_fun = PurpleAndYellow()
  ph <- pheatmap::pheatmap(sele.zscores, cluster_cols = F, 
                           cluster_rows = cluster.rows, 
                           show_colnames = F, fontsize = 10,
                           annotation_col = ann_column, 
                           color = color_fun,
                           annotation_colors = ann_colors, 
                           fontsize_row = 12)
  return(ph)
  
  
  
  
}


## virids style, downsample, match_cell, column color by sample & cluster
plot_enrich_tf2 <- function(sele.tfs, zscore.all, bc_clusters,
                            cluster.levels = NULL, sample.levels = NULL,
                            up.qt = 0.95, low.qt = 0.05,
                            ndownsample = 2000, 
                            match_cell = F, reScale = F, 
                            cluster.rows = F,
                            color_style = 'virid',
                            order.within = 'sample'){
  sele.zscores = zscore.all[sele.tfs, ]
  
  if(reScale) sele.zscores = t(scale(t(sele.zscores), center = T, scale = T))
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample and match_cell
  ncell.cl = min(table(bc_clusters$sample))
  set.seed(2020)
  bc_clusters.down = bc_clusters
  if(match_cell){
    bc_clusters.down = NULL
    for(sample0 in unique(bc_clusters$sample)){
      tmp = bc_clusters[sample == sample0]
      if(nrow(tmp) > ncell.cl) tmp = tmp[sort(sample((1:nrow(tmp)), ncell.cl)), ]
      bc_clusters.down = rbind(bc_clusters.down, tmp)
    }
  }
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = bc_clusters$cluster,
                          'barcode' = bc_clusters$barcode,
                          'sample' = bc_clusters$sample,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  nc = length(unique(bc_clusters$cluster))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(nc >= 3) color_cluster = getPalette(nc)
  if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
  names(color_cluster) = sort(unique(bc_clusters$cluster))
  
  nsample = length(unique(bc_clusters$sample))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(nsample <= 3) color_sample = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nsample]
  if(nsample > 3) color_sample = getPalette(nsample)
  names(color_sample) = sample.levels
  
  ## order sample by age
  ann_column = ann_column[order(factor(ann_column$sample,
                                       levels = sample.levels[sample.levels %in% ann_column$sample])), ]
  ann_column = ann_column[order(factor(ann_column$cluster,
                                       levels = sort(unique(ann_column$cluster)))), ]
  
  if(order.within == 'sample'){
    if(!is.null(cluster.levels)) {
      ann_column = ann_column[order(factor(ann_column$cluster,
                                           levels = cluster.levels)), ]
      
      color_cluster = color_cluster[cluster.levels]
    }
    
    if(!is.null(sample.levels)){
      ann_column = ann_column[order(factor(ann_column$sample,
                                           levels = sample.levels)), ]
      
      color_sample = color_sample[sample.levels]
    }
    
  }else{
    if(!is.null(sample.levels)){
      ann_column = ann_column[order(factor(ann_column$sample,
                                           levels = sample.levels)), ]
      
      color_sample = color_sample[sample.levels]
    }
    
    if(!is.null(cluster.levels)) {
      ann_column = ann_column[order(factor(ann_column$cluster,
                                           levels = cluster.levels)), ]
      
      color_cluster = color_cluster[cluster.levels]
    }
    
    
    
  }
  
  ann_colors = list('cluster' = color_cluster,
                    'sample' = color_sample)
  
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode <- NULL
  
  color_fun = viridis(100)
  if(color_style == 'purple-yellow') color_fun = PurpleAndYellow()
  
  ph <- pheatmap::pheatmap(sele.zscores, cluster_cols = F, 
                           cluster_rows = cluster.rows, 
                           show_colnames = F, fontsize = 10,
                           annotation_col = ann_column, 
                           color = color_fun,
                           annotation_colors = ann_colors, 
                           fontsize_row = 12)
  return(ph)
  
  
  
  
}

## prepare gene expression fraction in the corresponding cell populations ####
seurat.rna <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/42_add_finalized_cMeta_to_X01_samples/t.all.40.rds')

Ctypes = seurat.rna$ETP

mtx = seurat.rna@assays$RNA@counts

efreq.etp <- rowMeans(mtx[, Ctypes == 'ETP'] > 0)
efreq.non.etp <- rowMeans(mtx[, Ctypes == 'Non-ETP'] > 0)
efreq.near.etp <- rowMeans(mtx[, Ctypes == 'Near-ETP'] > 0)

save(efreq.etp, efreq.non.etp, efreq.near.etp, 
     file = 'MetaData2/expression_fraction_ETP.RData')



## prepare enriched TF list ####
chrvar.all <- readRDS('chromVARObj/chromVAR_40sample.rds')
zscore.all = chrvar.all@assays@data$z
dscore.all = chrvar.all@assays@data$deviations

## change to readable name
tfs = chrvar.all@elementMetadata$name
names(tfs) = NULL

rownames(dscore.all) = rownames(zscore.all) = tfs

save(zscore.all, dscore.all, 
     file = 'chromVARObj/zscore_dscore_allPeaks_ETP.RData')

load('chromVARObj/zscore_dscore_allPeaks_ETP.RData')
seurat.atac <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds')

ann.atac = data.table(seurat.atac@meta.data, keep.rownames = T)
cls = rep(c('ETP', 'Near-ETP', 'Non-ETP'),
          c(nrow(ann.atac[ETP == 'ETP']), nrow(ann.atac[ETP == 'Near-ETP']),
            nrow(ann.atac[ETP == 'Non-ETP'])))

cls = data.table('barcode' = c(ann.atac[ETP == 'ETP']$rn, 
                               ann.atac[ETP == 'Near-ETP']$rn, 
                               ann.atac[ETP == 'Non-ETP']$rn),
                 'cluster' = cls)

res.diff1 <- runDiffMotifEnrich(dscore.all, clusters = cls, fdr = 0.01,
                                max_cell_per_clust = 1000,
                                topn = 200, min_frac_per_cluster = 0.1)
res.diff1 = res.diff1[cluster1 == 'ETP']

res.diff2 <- runDiffMotifEnrich(dscore.all, clusters = cls, fdr = 0.01,
                                max_cell_per_clust = 1000,
                                topn = 200, min_frac_per_cluster = 0.1, control = 'ETP',
                                test.mode = 'control')
res.diff2 = res.diff2[cluster1 != 'ETP']

res.diff = rbind(res.diff1, res.diff2)
res.diff[, 'delta' := mean1 - mean0]


## filter by expression -- to do
res.diff = res.diff[delta > 0.01]
res.diff = res.diff[order(-delta)]

res.diff1 = res.diff[cluster1 == 'ETP' & feature %in% names(which(efreq.etp >= 0.1))]
res.diff2 = res.diff[cluster1 == 'Near-ETP' & feature %in% names(which(efreq.near.etp >= 0.1))]
res.diff3 = res.diff[cluster1 == 'Non-ETP' & feature %in% names(which(efreq.non.etp >= 0.1))]

res.diff = rbind(res.diff1, res.diff2, res.diff3)

ph.py <- plot_enrich_tf1(res.diff$feature, zscore.all, bc_clusters = cls,
                         up.qt = 0.97,
                         ndownsample = 1000, match_cell = T,
                         cluster.levels = c('Non-ETP', 'Near-ETP', 'ETP'), reScale = T,
                         color_style = 'purple-yellow')
ggsave(ph.py, filename = 'Figures/enriched_tfs_ETP_exprFracGT0.1.rds',
       width = 10, height = 10)

saveRDS(res.diff, file = 'MetaData2/enriched_tfs_ETP_exprFracGT0.1.rds')




## load DEGs list for the specific comparison ####
for(TRN.comparison in c('ETP_vs_Near-Non-ETP', 'Near-ETP_vs_ETP', 'Non-ETP_vs_ETP')){
  TRN.population = unlist(strsplit(TRN.comparison, '_'))[1]
  
  ## sele.type = all(default), stringent(using stringent DEGs and DA cutoffs), fixed (using selected TFs)
  sele.type = 'stringentPlus'  ## all, stringent or fixed or stringentPlus
  
  degs_list <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/16_construct-TRN/all.de.list.rds")
  degs = data.table(degs_list[[TRN.comparison]])
  degs = degs[abs(avg_log2FC) > 0.25 & p_val_adj < 0.01]
  
  if(sele.type == 'all') degs = degs[abs(avg_log2FC) > 0.5 & p_val_adj < 0.01]
  
  if(sele.type == 'stringent'){
    degs = degs[abs(avg_log2FC) > 0.5 & p_val_adj < 0.01]
    degs1 = degs[abs(avg_log2FC) > 1 & p_val_adj < 0.01 & group != 'de.tfs']
    degs2 = degs[abs(avg_log2FC) > 0.75 & p_val_adj < 0.01 & group == 'de.tfs']
    degs = rbind(degs1, degs2)
  }
  
  if(sele.type == 'stringentPlus'){
    degs = degs[abs(avg_log2FC) > 0.5 & p_val_adj < 0.01]
    degs1 = degs[abs(avg_log2FC) > 1 & p_val_adj < 0.01 & group != 'de.tfs']
    degs2 = degs[abs(avg_log2FC) > 0.75 & p_val_adj < 0.01 & group == 'de.tfs']
    
    surface.de = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/16_construct-TRN/de.surface.markers/de.surface.markers")
    dd.surface = data.table(surface.de[surface.de$group == 'de.genes', ])
    dd.surface = dd.surface[enriched ==  TRN.population, ]
    dd.surface[, c('comparison', 'enriched') := NULL]
    degs = rbind(degs1, degs2, dd.surface)
  }
  
  
  if(sele.type == 'fixed'){
    
    top.tfs.de.and.da = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/16_construct-TRN/DA_and_DE_TFs/top.tfs.de.and.da.rds")
    tfs = top.tfs.de.and.da[top.tfs.de.and.da$comparison == TRN.comparison, ]$gene
    degs1 = degs[abs(avg_log2FC) > 1 & p_val_adj < 0.01 & group != 'de.tfs']
    degs2 = degs[gene %in% tfs]
    degs = rbind(degs1, degs2)
  }
  
  ## construct TRN for a specify condition ####
  motif_ix = readRDS('MetaData2/motif_match_pv5E-05.rds')
  motif_name = motif_ix@colData$name
  pk_inf = motif_ix@rowRanges
  motif_pk_match <- motif_ix@assays@data$motifMatches
  rownames(motif_pk_match) = paste0(pk_inf@seqnames, '-', pk_inf@ranges)
  names(motif_name) = NULL
  colnames(motif_pk_match) = motif_name      
  
  
  ## < filter by degs ####
  predicted.ep = fread('EP_Prediction2/regrRes4_EP_overall.txt')
  
  predicted.ep = predicted.ep[gene_name %in% degs$gene]
  
  
  ## < filter by TF hits at enhancer side ####
  
  enriched.tf.dt = readRDS('MetaData2/enriched_tfs_ETP_exprFracGT0.1.rds')
  enriched.tf.dt = enriched.tf.dt[cluster1 == TRN.population, ]
  ntf = nrow(enriched.tf.dt)
  
  enriched.tfs = enriched.tf.dt$feature
  if(sele.type == 'stringent' || sele.type == 'stringentPlus'){
    enriched.tfs = enriched.tf.dt$feature[1:(min(10, ntf))]
  }
  if(sele.type == 'fixed'){
    enriched.tfs = NULL
  }
  
  # add tfs in degs
  degs.tf = degs[group == 'de.tfs' & avg_log2FC > 0]$gene
  enriched.tfs = unique(c(enriched.tfs, degs.tf))
  
  enriched.tf.dt = enriched.tf.dt[feature %in% enriched.tfs, ]
  enriched.tfs = intersect(enriched.tfs, colnames(motif_pk_match))
  
  enriched.tfs = enriched.tfs[!is.na(enriched.tfs)]
  sele.tf.mat = motif_pk_match[, enriched.tfs] 
  sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
  
  predicted.ep = predicted.ep[enhancer_peak %in% sele.peaks]
  
  
  
  if(T){
    ## filter by accessible peaks
    if(is.null(seurat.atac))seurat.atac <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds')
    
    ctypes = seurat.atac$ETP
    mtx = seurat.atac@assays$ATAC@counts
    frac.per.group <- sapply(c('ETP', 'Near-ETP',  'Non-ETP'), function(x){
      rowMeans(mtx[, ctypes == x] > 0)
    })
    frac.max.peak <- frac.per.group[, TRN.population]
    
    final.peaks.str = names(which(frac.max.peak > 0.1))
    final.peaks.str = sapply(final.peaks.str, function(x) unlist(strsplit(x, ','))[1])
    names(final.peaks.str) = NULL
    predicted.ep = predicted.ep[enhancer_peak %in% final.peaks.str ]
  }
  
  
  ## < construct network ####
  ## split loop by tf
  ep.tf = list()
  for(TF0 in enriched.tfs){
    peaks0 = names(which(sele.tf.mat[, TF0] > 0))
    ep0 = predicted.ep[enhancer_peak %in% peaks0]
    ep0 = subset(ep0, select = c('gene_name', 'Estimate', 'fdr'))
    ep0[, 'score' := -log10(fdr)]
    ep0$TF = TF0
    ep0 = ep0[score > 2]
    ep0 = rbind(ep0[gene_name %in% degs[avg_log2FC > 0]$gene & Estimate > 0], 
                ep0[gene_name %in% degs[avg_log2FC < 0]$gene & Estimate < 0])
    ep0[, 'N' := .N, by = gene_name]
    ep0[, 'score' := sum(score), by = gene_name]
    ep0 = subset(ep0, select = c('TF', 'gene_name', 'score'))
    ep.tf[[TF0]] = ep0[!duplicated(ep0)]
  }
  
  ep.tf.comb = do.call('rbind', ep.tf)
  
  ep.tf.comb = ep.tf.comb[order(-score)]
  #ep.tf.comb = ep.tf.comb[1:100, ]
  
  vertex.gr = data.table('gene_name' = unique(c(ep.tf.comb$TF, ep.tf.comb$gene_name)),
                         'group' = 'Gene')
  vertex.gr$group[vertex.gr$gene_name %in% enriched.tfs] = 'TF'
  setkey(vertex.gr, gene_name)
  
  grn <- graph_from_edgelist(as.matrix(ep.tf.comb[, 1:2]), directed = F)
  
  V(grn)$vgroup = vertex.gr[J(V(grn)$name)]$group 
  
  
  plot.igraph(grn,  vertex.color=c( "skyblue", "red")[1+(V(grn)$vgroup=="TF")],
              vertex.label.dist=0, layout = layout_with_graphopt,
              vertex.shape = c('square', 'circle')[1+(V(grn)$vgroup=="TF")],
              margin = 0)
  
  
  ## add expression information
  degs0 = copy(degs)
  degs0[, 'group' := NULL]
  degs0 %<>% unique()
  
  vertex.gr[gene_name %in% degs0$gene, 
            'expr_logFC' := degs0[gene == gene_name]$avg_log2FC, by = gene_name]
  
  vertex.gr[gene_name %in% degs0$gene, 
            'NegLog10PV' := -log10(degs0[gene == gene_name]$p_val_adj), by = gene_name]
  
  vertex.gr[gene_name %in% enriched.tf.dt$feature, 
            'tf_dev' := abs(enriched.tf.dt[feature == gene_name]$delta), by = gene_name]
  vertex.gr[gene_name %in% enriched.tfs, 
            'NegLog10PV' := -log10(enriched.tf.dt[feature == gene_name]$pv_adjust), by = gene_name]
  
  vertex.gr$effect_size = vertex.gr$tf_dev /max(vertex.gr$tf_dev, na.rm = T)
  vertex.gr$tmp = vertex.gr$expr_logFC /max(abs(vertex.gr$expr_logFC), na.rm = T)
  vertex.gr[, 'effect_size' := ifelse(is.na(tf_dev), tmp, effect_size), by = gene_name]
  vertex.gr[, 'tmp' := NULL]
  vertex.gr[is.infinite(NegLog10PV)]$NegLog10PV = max(vertex.gr[!is.infinite(NegLog10PV)]$NegLog10PV)+10
  
  ep.tf.comb$direction = 'up'
  ep.tf.comb[gene_name %in% degs0[avg_log2FC < 0]$gene]$direction = 'down'
  
  vertex.gr$direction = 'up'
  
  vertex.gr[gene_name %in% degs0[avg_log2FC < 0]$gene]$direction = 'down'
  
  
  vertex.gr[, 'score' := ifelse(direction == 'down', -NegLog10PV, NegLog10PV)]
  
  vfile = paste0('EP_Prediction2/TRN_nodes_enrichedTFs_in_', TRN.population, '.txt')
  efile = paste0('EP_Prediction2/TRN_edges_enrichedTFs_in_', TRN.population, '.txt')
  
  if(sele.type != 'all'){
    vfile = paste0('EP_Prediction2/TRN_nodes_enrichedTFs_in_', TRN.population, '_', sele.type, '.txt')
    efile = paste0('EP_Prediction2/TRN_edges_enrichedTFs_in_', TRN.population, '_', sele.type, '.txt')
  }
  
  write.table(ep.tf.comb, file = efile,
              sep = '\t', row.names = F, quote = F)
  write.table(vertex.gr, file = vfile,
              sep = '\t', row.names = F, quote = F)
  
}
