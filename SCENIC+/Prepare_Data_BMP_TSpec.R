library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(anndata)
library(Signac)
library(Matrix)
library(ggplot2)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL/SCENICPlus")

rna = readRDS("etp.25.g1.downsample.1711.each.patient.rds")
table(rna$comparison.bmp.vs.t.specified)
table(rna$ETP)

#Second for ETP 
rna = readRDS("../GEO_Upload/Seurat/Full_CITE_seq_40_TALL.rds")
table(rna$ETP)
table(rna$orig.ident)
table(rna$is.blast.viscello)
rna = subset(rna, subset = is.blast.viscello == TRUE)
Idents(rna) = "orig.ident"
rna.subtype = subset(rna, downsample = 1146)
table(rna.subtype$ETP)
table(rna.subtype$orig.ident)
table(rna.subtype$is.blast.viscello)

atac.blasts = readRDS("t.all.40.atac.blasts.1146.each.with.v4.anno.updated.RDS")
table(atac.blasts$comparison.bmp.vs.t.specified)
table(atac.blasts$ETP)

atac.blasts = subset(atac.blasts, comparison.bmp.vs.t.specified == "Other", invert = T)
table(atac.blasts$comparison.bmp.vs.t.specified)
table(atac.blasts$ETP)

DimPlot(atac.blasts, group.by = "comparison.bmp.vs.t.specified") + coord_fixed()

########
#Save data for SCENIC+
########
#Save RNA-seq
DefaultAssay(rna) <- "RNA"
rna.loom <- as.loom(rna, filename = "Data_SCENICplus/TALL_BMP_TSpec_RNA.loom", verbose = TRUE) 
rna.loom$close_all() 

DefaultAssay(rna.subtype) <- "RNA"
rna.loom <- as.loom(rna.subtype, filename = "Data_SCENICplus/TALL_ETP_Subtype_RNA.loom", verbose = TRUE) 
rna.loom$close_all() 

#Save ATAC-seq data 
#Sparse matrix 
counts_matrix <- as.matrix(atac.blasts@assays$ATAC@counts)
counts_sparse <- Matrix::Matrix(counts_matrix , sparse = T)
writeMM(counts_sparse, file = "Data_SCENICplus/ATAC_Peaks_Sparse.mtx")
cell_names <- colnames(counts_sparse)
region_names <- rownames(counts_sparse)
cell_names_file <- "Data_SCENICplus/ATAC_Cell_Names.txt"
region_names_file <- "Data_SCENICplus/ATAC_Region_Names.txt"
metadata_frame <- atac.blasts@meta.data
write.table(cell_names, file = cell_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(region_names, file = region_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(metadata_frame, file = "Data_SCENICplus/ATAC_Metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Save regions to bed format 
convert_to_bed <- function(region) {
  parts <- unlist(strsplit(region, "-"))
  chr <- parts[1]
  start <- as.integer(parts[2])
  end <- as.integer(parts[3])
  bed <- c(chr, start, end)
  return(bed)
}
bed_data <- apply(as.matrix(region_names), 1, convert_to_bed)
write.table(t(bed_data), "Data_SCENICplus/ATAC_Region_Names.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
