## using all normal cell types as reference
## test on blast cells

library(infercnv)
library(Seurat)
library(data.table)


seurat.T.myeloid = readRDS('/mnt/isilon/tan_lab/chenc6/ETP_ALL/CITE-Seq/Scripts/SeuratObject/2_9_blast_healthy_anno_10samples_ETPrelapsed.rds')
Tmyeloid_samples = unique(seurat.T.myeloid$sample.name)

# load annotation data
Tmyeloid_ann = seurat.T.myeloid@meta.data

seurat.obj <- subset(seurat.T.myeloid, sample.name == "PATEVG")
Tmyeloid_ann <- subset(Tmyeloid_ann, sample.name == "PATEVG")
rm(seurat.T.myeloid)

if(any(rownames(Tmyeloid_ann) != colnames(seurat.obj))) stop('Cell name not matched between seurat obj and annotation!')
mtx <- seurat.obj@assays$RNA@counts
rm(seurat.obj)
Tmyeloid_ann = subset(Tmyeloid_ann, select = level.1.anno)

Tmyeloid_ann$level.1.anno = ifelse(grepl(Tmyeloid_ann$level.1.anno, pattern = 'blast'),
                        'blast', Tmyeloid_ann$level.1.anno)
Tmyeloid_ann$level.1.anno = ifelse(!grepl(Tmyeloid_ann$level.1.anno, pattern = 'blast'),
                        'normal', Tmyeloid_ann$level.1.anno)
mtx = mtx[, Tmyeloid_ann$level.1.anno != 'NA']
Tmyeloid_ann = subset(Tmyeloid_ann, level.1.anno != 'NA')

exclude_minors_ctypes = names(which(table(Tmyeloid_ann$level.1.anno) < 10))
mtx = mtx[, !Tmyeloid_ann$level.1.anno %in% exclude_minors_ctypes]
Tmyeloid_ann = subset(Tmyeloid_ann, !level.1.anno %in% exclude_minors_ctypes)

names(Tmyeloid_ann) <- NULL


genelist = read.table('/mnt/isilon/tan_lab/chenc6/Tools/InferCNV/Hg38_Gene_Position.txt', 
                      sep = '\t', stringsAsFactors = F)

genelist = genelist[!duplicated(genelist$V1), ]

rownames(genelist) <- genelist$V1

#reorder chr
chrs = paste0('chr', 1:22)
genes_ann = NULL
for(chr0 in chrs) genes_ann = rbind(genes_ann, genelist[genelist$V2 == chr0, ])

genes_ann$V1 <- NULL

names(genes_ann) <- NULL



infercnv_obj = CreateInfercnvObject(raw_counts_matrix = as.matrix(mtx),
                                    annotations_file = Tmyeloid_ann,
                                    delim = "\t",
                                    gene_order_file = genes_ann,
                                    ref_group_names = 'normal')

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             min_cells_per_gene = 5,
                             out_dir = paste0('/mnt/isilon/tan_lab/chenc6/ETP_ALL/CITE-Seq/Scripts/InferCNV_Objects_blast/PATEVG_new1'), 
                             cluster_by_groups = FALSE, 
                             k_obs_groups = 4,
                             analysis_mode = "subclusters",
                             denoise = T,
                             HMM = TRUE)

saveRDS(infercnv_obj, file = paste0('/mnt/isilon/tan_lab/chenc6/ETP_ALL/CITE-Seq/Scripts/InferCNV_Objects_blast/PATEVG_new1/final_obj_PATEVG_new1.rds'))



infercnv.obj <- readRDS("InferCNV_Objects_blast/PAVVVK/final_obj_PAVVVK.rds")
obs_hcl <- infercnv.obj@tumor_subclusters$hc[["blast"]]
split_groups <- cutree(obs_hcl, k=2)

split_groups <- data.frame(split_groups)

seurat.obj <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/other_r_objects/combined_seurat_objects/2021-1-21_34samples_with_nn_annos_ref_umap.rds")

seurat.obj.PAVVVK <- subset(seurat.obj, sample.name == "PAVVVK")

seurat.obj.PAVVVK <- AddMetaData(seurat.obj.PAVVVK, split_groups)

DimPlot(seurat.obj.PAVVVK, reduction = 'umap_50', group.by = "split_groups")

infercnv.obj@tumor_subclusters = NULL
plot_cnv(
  infercnv.obj,
  out_dir = ".",
  title = "inferCNV",
  obs_title = "Observations (Cells)",
  ref_title = "References (Cells)",
  cluster_by_groups = FALSE,
  cluster_references = TRUE,
  k_obs_groups = 3,
  contig_cex = 1,
  x.center = mean(infercnv_obj@expr.data),
  x.range = "auto",
  hclust_method = "ward.D",
  color_safe_pal = FALSE,
  output_filename = "infercnv",
  output_format = "png",
  png_res = 300,
  dynamic_resize = 0,
  ref_contig = NULL,
  write_expr_matrix = FALSE
)


infercnv.obj <- readRDS("InferCNV_Objects_blast/Replot/PAURXZ/final_obj_PAURXZ.rds")
obs_hcl <- infercnv.obj@tumor_subclusters$hc[["blast"]]
split_groups <- cutree(obs_hcl, k=2)

group <- fread("InferCNV_Objects_blast/Replot/PAURXZ/infercnv_PAURXZ.observation_groupings.txt")

row.names(group) <- group$V1

split_groups <- data.frame(split_groups)

seurat.obj <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/other_r_objects/combined_seurat_objects/2021-1-21_34samples_with_nn_annos_ref_umap.rds")

seurat.obj.PAURXZ <- subset(seurat.obj, sample.name == "PAURXZ")

seurat.obj.PAURXZ <- AddMetaData(seurat.obj.PAURXZ, group)

png("InferCNV_Objects_blast/Replot/PAURXZ/UMAP_ref_PAURXZ.png", width = 900, height = 600)
DimPlot(seurat.obj.PAURXZ, reduction = 'umap_50', group.by = "Dendrogram Group")
dev.off()
