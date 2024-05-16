library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(GenomicRanges)

## check motif overlapping with any peak ####
#seurat.atac = readRDS('SeuratObj/Seurat_25ETPsample_downsample_VFACS.rds')
seurat.atac <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds')

peaks = rownames(seurat.atac)
peaks = sapply(peaks, function(x) unlist(strsplit(x, ','))[1])
names(peaks) = NULL

pks = tidyr::separate(data = data.table('peak_name' = peaks),
                      col = 'peak_name', into = c('chr', 'start', 'end'),
                      remove = F)
pks$start = as.integer(pks$start)
pks$end = as.integer(pks$end)
setkey(pks, chr, start)
# Make a set of peaks
peaks <- GenomicRanges::GRanges(seqnames = pks$chr,
                 ranges = IRanges::IRanges(start = pks$start,
                                  end = pks$end))

motif_ix <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = 1e-3)
motif_ix_str <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = 5e-5)
if(F){
        motif_name = motif_ix@colData$name
        pk_inf = motif_ix@rowRanges
        motif_pk_match <- motif_ix@assays@data$motifMatches
        rownames(motif_pk_match) = paste0(pk_inf@seqnames, '-', pk_inf@ranges)
        names(motif_name) = NULL
        colnames(motif_pk_match) = motif_name      
}

#saveRDS(motif_ix, file = 'MetaData/motif_match_pvE-03.rds')
#saveRDS(motif_ix_str, file = 'MetaData/motif_match_pv5E-05.rds')

saveRDS(motif_ix, file = 'MetaData2/motif_match_pvE-03.rds')
saveRDS(motif_ix_str, file = 'MetaData2/motif_match_pv5E-05.rds')
