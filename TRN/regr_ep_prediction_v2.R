source('scDataAnalysis_Utilities_simp.R')
`%notin%` = Negate(`%in%`)

## get a binary matrix indicates the gene-peak affinity
## gene.list are data.table including column name gene_name 
## gene_ann should include gene_name,chr,start,end
get_gene2peak_map <- function(gene.list, peak_names,
                              gene_ann, distal_dist = 2e05){
  
  peaks <- tidyr::separate(data.table(peak_name = peak_names), 
                           col = peak_name, remove = F,
                           into = c('chr', 'start', 'end'))
  class(peaks$start) = class(peaks$end) = 'integer'
  setkey(gene_ann, gene_name)
  setkey(peaks, peak_name)
  gene.list = gene.list[gene_name %in% gene_ann$gene_name]
  gene.list$chr = gene_ann[gene.list$gene_name]$chr
  gene.list$start = gene_ann[gene.list$gene_name]$start
  gene.list$end = gene_ann[gene.list$gene_name]$end
  
  ## for each gene, get the corresponding peaks
  gene2peaks = lapply(gene.list$gene_name, function(x) {
    
    chr0 = gene_ann[x]$chr
    start0 = gene_ann[x]$start
    end0 = gene_ann[x]$end
    
    peaks0 = peaks[chr == chr0]
    peaks0 = peaks0[abs(start/2 + end/2 - start0/2 - end0/2) <= distal_dist]
    return(peaks0$peak_name)
  } )
  
  ## pool all peaks relate to one gene
  gene2peaks.u <- lapply(sort(unique(gene.list$gene_name)), function(x){
    id = which(gene.list$gene_name == x)
    tmp_peak <- do.call('c', lapply(id, function(x) gene2peaks[[x]]))
    return(tmp_peak)
  })
  names(gene2peaks.u) <- sort(unique(gene.list$gene_name))
  lens = sapply(gene2peaks.u, length)
  
  genes.f <- names(which(lens > 0))
  lens = lens[lens > 0]
  ## construct overlap matrix
  gene2peaks.dt <- data.table('gene' = rep(genes.f, lens),
                              'peak' = do.call('c', lapply(genes.f, 
                                                           function(x) gene2peaks.u[[x]])))
  upeaks = sort(unique(gene2peaks.dt$peak))
  gene2peaks.dt[, 'id1' := which(genes.f == gene), by = gene]
  gene2peaks.dt[, 'id2' := which(upeaks == peak), by = peak]
  gene2peak.map <- sparseMatrix(i = gene2peaks.dt$id1,
                                j = gene2peaks.dt$id2,
                                dimnames = list(genes.f, upeaks))
  gene2peak.map = gene2peak.map * 1
  
  return(gene2peak.map)
}

## annotate peaks with gene +/- 5kb of its TSS
# input peak_coords with chr-start-end, format
annPeak2Gene <- function(peak_coords, gene_ann, proxim_dist = 5e+03){
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  peaks = tidyr::separate(data.table(x = peak_coords),
                          col = x,
                          into = c('chr', 'start', 'end'))
  
  
  peaks$peak_name = peak_coords
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  chrs = unique(peaks$chr)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    genes0 = gene_ann[chr == chr0]
    
    peaks0$gene_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= (peaks0$end[i] + proxim_dist) & tss >= (peaks0$start[i] - proxim_dist)]
      if(nrow(tss0) > 0 ) {
        peaks0$gene_name[i] = paste(unique(tss0$gene_name), collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  setkey(peaks_ann, peak_name)
  
  return(peaks_ann)
  
}


## map gene to overlapping atac peak
## gene_list with genename, chr, start, end
geneOverlapPeak <- function(gene_list, peak_names, mid_dist = 1000){
  # should include tss information in gene_list
  peaks = tidyr::separate(data = data.table('peak_name' = peak_names),
                          col = peak_name, into = c('chr', 'start', 'end'),
                          remove = F)
  class(peaks$chr) = 'character'
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  chrs = unique(gene_list$chr)
  gene_new = NULL
  peaks[, 'midP' := start/2 + end/2]
  for(chr0 in chrs){
    gene0 = gene_list[chr == chr0, ]
    gene0$peak_name = 'Not_Found'
    peaks0 = peaks[chr == chr0]
    gene0[, 'peak_id0' := any( abs(peaks0$midP -start) < mid_dist | abs(peaks0$midP - end) < mid_dist), 
          by = gene_name]
    gene1 = gene0[peak_id0 == FALSE]
    gene2 = gene0[peak_id0 == TRUE]
    gene2[, 'peak_id' :=  which.min(abs(peaks0$midP - start - 1000)), by = gene_name]
    
    gene2[, 'peak_name' :=  peaks0[peak_id]$peak_name, by = gene_name]
    gene2$peak_id = NULL
    gene_new = rbind(gene_new, gene1, gene2)
  }
  gene_new[, c('peak_id0') := NULL]
  return(gene_new)
}


## combine all coembed results from a single patient ####
seurat.atac <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds')

sampleNames = unique(seurat.atac$sample.name)

final_matching = smooth.rna = smooth.atac = list()
for(sampleName in sampleNames){
  final_matching[[sampleName]] <- readRDS(paste0('Coembed_Results2/atac_rna_coembedding_cell_matching_',
                                                 sampleName, '.rds'))
  smooth.rna[[sampleName]] <- as.data.frame(t(readRDS(paste0('Coembed_Results2/rna_metacell_expr_',
                                                 sampleName, '.rds'))))
  smooth.atac[[sampleName]] <- as.data.frame(t(readRDS(paste0('Coembed_Results2/atac_metacell_access_',
                                             sampleName, '.rds'))))
}

final_matching = do.call('rbind', final_matching)
smooth.rna = data.table::rbindlist(smooth.rna)
smooth.atac = data.table::rbindlist(smooth.atac)
smooth.rna$cell_name = final_matching$rna_cell
smooth.atac$cell_name = final_matching$atac_cell

saveRDS(final_matching, file = 'Coembed_Results2/combined_final_matching.rds')
saveRDS(smooth.rna, file = 'Coembed_Results2/combined_scRNA_smoothed.rds')
saveRDS(smooth.atac, file = 'Coembed_Results2/combined_scATAC_smoothed.rds')


## regression ####
final_matching = readRDS('Coembed_Results2/combined_final_matching.rds')
smooth.rna = readRDS('Coembed_Results2/combined_scRNA_smoothed.rds')
smooth.atac = readRDS('Coembed_Results2/combined_scATAC_smoothed.rds')

## filter peaks that accessible in less than 5% of all cell type
if(T){
  
  peaks.mean.ctype <- sapply(unique(seurat.atac$ETP), function(x){
    rowMeans(seurat.atac@assays$ATAC@data[, seurat.atac$ETP == x] > 0)
  })
  
  rmax = apply(peaks.mean.ctype, 1,  max)
  summary(rmax)
  filtered.peaks = names(which(rmax > 0.05))
  filtered.peaks = lapply(filtered.peaks, function(x) unlist(strsplit(x, ','))[1])
  filtered.peaks = do.call('c', filtered.peaks)
}

cand.peaks = colnames(smooth.atac)
cand.peaks = sapply(cand.peaks, function(x) unlist(strsplit(x, ','))[1])
names(cand.peaks) = NULL
colnames(smooth.atac) = cand.peaks
all(filtered.peaks %in% cand.peaks)


atac.cnames = smooth.atac$cell_name
smooth.atac = subset(smooth.atac, select = filtered.peaks)
acces.frac = sapply(filtered.peaks, function(x) mean(smooth.atac[[x]] > 0))

final.peaks = filtered.peaks[acces.frac > 0.01]
smooth.atac = subset(smooth.atac, select = final.peaks)


## focus on selected genes (or degs)
degs_list = all.de.list <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/16_construct-TRN/all.de.list.rds")
degs = unique(c(degs_list[[1]]$gene, degs_list[[2]]$gene, degs_list[[3]]$gene))
surface.de = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/16_construct-TRN/de.surface.markers/de.surface.markers")
degs = unique(c(degs), unique(surface.de[surface.de$group == 'de.genes',]$gene))
## gene-peak-affinity map
## construct gene peak affinity binary matrix
gene_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/gene_ann_hg38.txt')
gene_ann[, 'Tss' := ifelse(strand == '+', start, end)]

gene2peak.map <- get_gene2peak_map(gene.list = data.table('gene_name' = degs),
                                   peak_names = final.peaks, 
                                   gene_ann = gene_ann,
                                   distal_dist = 5e5)

peaks.used = colnames(gene2peak.map)
degs = degs[degs %in% rownames(gene2peak.map)]

## do regression gene by gene
smooth.rna = subset(smooth.rna, select = degs)
smooth.atac = subset(smooth.atac, select = peaks.used)

stime = Sys.time()
regr.list = list()
for(x in degs){
  exprs <- as.numeric(smooth.rna[[x]])
  names(exprs) <- NULL
  peaks = names(which(gene2peak.map[x, ] == 1))
  covrs = as.matrix(subset(smooth.atac, select = peaks))
  
  rdata <- data.frame(cbind(exprs, covrs))
  res <- coef(summary(lm(exprs ~ ., data = rdata)))
  colnames(res)[4] <- 'P_value'
  regr.list[[x]] <- res
}
etime = Sys.time()
etime-stime
names(regr.list) = degs
saveRDS(regr.list, "EP_Prediction2/regrRes4ep_prediction.rds")



## summarize/filter loops ####
regr.list = readRDS("EP_Prediction2/regrRes4ep_prediction.rds")
regr.sum <- lapply(degs, function(t){
  x = regr.list[[t]]
  x = x[, c(1, 4)]
  x = data.frame(x)
  x = data.table(x, keep.rownames = T)
  x$gene_name = t
  return(x)
})

regr.sum = do.call('rbind', regr.sum)

regr.sum[, 'p_val_adj' := pmin(1, P_value*nrow(regr.sum))]
regr.sum$fdr = p.adjust(regr.sum$P_value, method = 'fdr')

regr.filtered = regr.sum[fdr < 0.05 & abs(Estimate) > 0.3 & grepl(rn, pattern = '^chr')]

regr.filtered$peak_name = sapply(regr.filtered$rn, function(x) gsub('.', '-', x, fixed = T)  )
regr.filtered$rn <- NULL

## to visualize on ucsc genome browser
## (promoter side: closest peak to TSS -- gene level)
gene_ann.deg = gene_ann[gene_name %in% regr.filtered$gene_name, ]

gene_ann.deg[, 'promoter_start' := Tss - 1000]
gene_ann.deg[, 'promoter_end' := Tss + 1000]

setkey(gene_ann.deg, gene_name)
regr.filtered[, 'start' := gene_ann.deg[J(regr.filtered$gene_name)]$promoter_start]
regr.filtered[, 'end' := gene_ann.deg[J(regr.filtered$gene_name)]$promoter_end]
regr.filtered[, 'chr' := gene_ann.deg[J(regr.filtered$gene_name)]$chr]
regr.filtered[, 'promoter_pos' := paste(chr, start, end, sep = '-')]

## filter otherend not overlapping with promoters
tss_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/transcript_ann_hg38.txt')
tss_ann = tss_ann[gene_biotype %in% c('protein_coding', 'lincRNA', 'miRNA')]
tss_ann[, 'Tss' := ifelse(strand == '+', start, end)]

# any peak within promoter region
peak.ann = annPeak2Gene(peaks.used, tss_ann, 2000)
setkey(peak.ann, peak_name)
peaks.nprom = peak.ann[nchar(gene_name) == 0]$peak_name
peaks.prom = peak.ann[nchar(gene_name) > 0]$peak_name

regr.filtered.ep = regr.filtered[peak_name %in% peaks.nprom]

## assign nearest peak to promoter
gene_list <- subset(regr.filtered.ep, 
                    select = c(gene_name, chr, start, end)) %>%
  .[!duplicated(.)]

gene_list2peak = geneOverlapPeak(gene_list, peak_names = peaks.prom,
                                 mid_dist = 1000)
gene_list2peak = gene_list2peak[peak_name != 'Not_Found']
setkey(gene_list2peak, gene_name)

regr.filtered.ep = regr.filtered.ep[gene_name %in% gene_list2peak$gene_name]
regr.filtered.ep[, 'promoter_peak' := gene_list2peak[J(regr.filtered.ep$gene_name)]$peak_name]


regr.filtered.ep = subset(regr.filtered.ep, select = c(gene_name, promoter_pos, 
                                                       promoter_peak, peak_name,
                                                       P_value, p_val_adj, fdr, Estimate))
names(regr.filtered.ep)[4] = 'enhancer_peak'

regr.filtered.ep[, 'chr1' := unlist(strsplit(promoter_peak, '-'))[1], by = promoter_peak]
regr.filtered.ep[, 'start1' := as.integer(unlist(strsplit(promoter_peak, '-'))[2]), by = promoter_peak]
regr.filtered.ep[, 'end1' := as.integer(unlist(strsplit(promoter_peak, '-'))[3]), by = promoter_peak]

regr.filtered.ep[, 'chr2' := unlist(strsplit(enhancer_peak, '-'))[1], by = enhancer_peak]
regr.filtered.ep[, 'start2' := as.integer(unlist(strsplit(enhancer_peak, '-'))[2]), by = enhancer_peak]
regr.filtered.ep[, 'end2' := as.integer(unlist(strsplit(enhancer_peak, '-'))[3]), by = enhancer_peak]

regr.filtered.ep[, 'ep_dist' := abs(start1 + end1 - start2 - end2)/2]
regr.filtered.ep = subset(regr.filtered.ep, select = c(gene_name, promoter_pos, 
                                                       promoter_peak, enhancer_peak, ep_dist,
                                                       P_value, p_val_adj, fdr, Estimate))

fwrite(regr.filtered.ep, file = 'EP_Prediction2/regrRes4_EP_overall.txt',
       sep = '\t')
