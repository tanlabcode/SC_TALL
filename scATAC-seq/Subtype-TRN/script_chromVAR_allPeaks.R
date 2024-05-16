library(Seurat)

##reture chromvar result given a seurat obj

source('scDataAnalysis_Utilities_simp.R')
args = commandArgs(T)
seuratPath = args[1]
chromvarPath = args[2]
ndownSample = as.integer(args[3])

seurat.obj = readRDS(seuratPath)

mtx = seurat.obj@assays$ATAC@counts

## further filter peaks
rs = Matrix::rowSums(mtx > 0)
mtx = mtx[rs > 10, ]

rn = rownames(mtx)
new_rn = sapply(rn, function(x) unlist(strsplit(x, ','))[1])

rownames(mtx) = new_rn
if(ndownSample < ncol(mtx)) {
    set.seed(2020)
    sele.cells = sample((1:ncol(mtx)), ndownSample)
    mtx = mtx[, sele.cells]
}
chromVAR.obj = run_chromVAR(mtx, genomeName = 'BSgenome.Hsapiens.UCSC.hg38')

saveRDS(chromVAR.obj, file = chromvarPath) 
