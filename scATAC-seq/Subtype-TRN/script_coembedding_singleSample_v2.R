source('scDataAnalysis_Utilities_simp.R')
library(ggplot2)
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
  idf = log(1 + idf)
  norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
  norm.data[which(x = is.na(x = norm.data))] <- 0
  return(norm.data)
}


seurat.rna <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/42_add_finalized_cMeta_to_X01_samples/t.all.40.rds')
seurat.atac <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds')

args = commandArgs(T)
sampleName = args[1]

seurat.rna0 = subset(seurat.rna, sample.name == sampleName)
seurat.atac0 = subset(seurat.atac, sample.name == sampleName)

## downsample to 4k if more --rna
if(ncol(seurat.rna0) > 4000){
  set.seed(2021)
  seurat.rna0$bc = colnames(seurat.rna0)
  sele.cells = sample(seurat.rna0$bc, 4000)
  seurat.rna0 = subset(seurat.rna0, bc %in% sele.cells)
  seurat.rna0$bc <- NULL
}
rm(seurat.rna, seurat.atac)

## do standard seurat pipeline, but with pearson residual normalization v2

DefaultAssay(seurat.rna0) <- 'RNA'
seurat.rna0 = FindVariableFeatures(seurat.rna0)
genes4anchors <- VariableFeatures(object = seurat.rna0)
seurat.rna0 <- ScaleData(seurat.rna0)
seurat.rna0 = RunPCA(seurat.rna0, npcs = 50, verbose = F)
seurat.rna0 = RunUMAP(seurat.rna0, reduction = 'pca', dims = 1:50)


atac.mtx = seurat.atac0@assays$ATAC@counts
seurat.atac0@assays$ATAC@data <- TF_IDF(atac.mtx, F)

rn = rownames(atac.mtx)
rownames(atac.mtx) <- sapply(rn, function(x) unlist(strsplit(x, ','))[1])
activity.matrix = generate_gene_cisActivity('/mnt/isilon/tan_lab/chenc6/Tools/SingleCellAnnotation/GRCh38/genes/genes.gtf',
                                            atac.mtx, 
                                            include_body = T)

seurat.atac0[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

DefaultAssay(seurat.atac0) <- "ACTIVITY"
seurat.atac0 <- NormalizeData(seurat.atac0)
seurat.atac0 <- FindVariableFeatures(seurat.atac0)
VariableFeatures(seurat.atac0) = genes4anchors
seurat.atac0 <- ScaleData(seurat.atac0)

DefaultAssay(seurat.atac0) <- "ATAC"
seurat.atac0$tech = 'ATAC'
seurat.rna0$tech = 'RNA'

## transfer label 

#genes4anchors = NULL
transfer.anchors <- FindTransferAnchors(reference = seurat.rna0,
                                        query = seurat.atac0,
                                        features = genes4anchors,
                                        reference.assay = "RNA",
                                        normalization.method = "LogNormalize",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5)

#co-embedding
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to

refdata <- GetAssayData(seurat.rna0, assay = "RNA", slot = "data")[genes4anchors, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat.atac0[["pca"]],
                           dims = 1:ncol(seurat.atac0[["pca"]]))

# this line adds the imputed data matrix to the seurat.atac0 object
seurat.atac0[["RNA"]] <- imputation
coembed <- merge(x = seurat.rna0, y = seurat.atac0)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes4anchors, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes4anchors, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:50)
p0 <- DimPlot(coembed, group.by = 'tech', label = T) 
ggsave(p0, file = paste0('Coembed_Results2/coembed_', sampleName, '.png'), device = 'png')
#saveRDS(coembed, file = 'Coembed_Results/seurat_50Kcoembed_30ETP.rds')

## 1 to 1 cell matching ####
umap_coproj = coembed@reductions$umap@cell.embeddings
ac_cells <- colnames(coembed)[coembed$tech == "ATAC"]
rna_cells <- colnames(coembed)[coembed$tech == "RNA"]
umap.rna = umap_coproj[rna_cells, ]
umap.atac = umap_coproj[ac_cells, ]
final_matching <- data.table(atac_cell = ac_cells)
final_matching$atac_cell <- as.character(final_matching$atac_cell)

dist0 <- pracma::distmat(umap.atac, umap.rna)
final_matching$rna_cell <- sapply(1:nrow(umap.atac), function(x) 
  names(which.min(dist0[x, ])))

saveRDS(final_matching, paste0("Coembed_Results2/atac_rna_coembedding_cell_matching_", sampleName,
                               ".rds"))


## smoothing ####

K = 10
seurat.rna0 <- FindNeighbors(seurat.rna0, k.param = K, reduction = 'pca')
knn.mat.rna = (seurat.rna0@graphs$RNA_nn > 0) * 1

seurat.atac0 <- FindNeighbors(seurat.atac0, k.param = K, reduction = 'pca')
knn.mat.atac = (seurat.atac0@graphs$ATAC_nn > 0) * 1

all(rowSums(knn.mat.rna) == K)
all(rowSums(knn.mat.atac) == K)

smooth.rna = seurat.rna0@assays$RNA@data %*% t(knn.mat.rna)
smooth.atac = seurat.atac0@assays$ATAC@data %*% t(knn.mat.atac)

smooth.rna = smooth.rna[, final_matching$rna_cell]
smooth.atac = smooth.atac[, final_matching$atac_cell]

saveRDS(smooth.rna, file = paste0("Coembed_Results2/rna_metacell_expr_", sampleName, ".rds"))
saveRDS(smooth.atac, file = paste0("Coembed_Results2/atac_metacell_access_", sampleName, ".rds"))

