library(biomaRt)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
library(cicero)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(SeuratWrappers)
library(JASPAR2022)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RColorBrewer)
library(future)
library(gridExtra)
library(ggpubr)
library(speckle)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(monocle)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL")

atac.blasts = readRDS("t.all.40.atac.blasts.1146.each.with.v4.anno.rds")

atac_counts = GetAssayData(atac.blasts, assay = "ATAC", slot = "counts")
atac_data = GetAssayData(atac.blasts, assay = "ATAC", slot = "data")
features <- sub(",.*", "", rownames(atac.blasts))
rownames(atac_counts) = features
rownames(atac_data) = features

atac.blasts[["ATAC"]] = CreateChromatinAssay(counts = atac_counts, genome = "hg38")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(atac.blasts) <- annotations
#saveRDS(atac.blasts, file = "t.all.40.atac.blasts.1146.each.with.v4.anno.updated.RDS")

atac.cds <- as.CellDataSet(x = atac.blasts)
atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDimS(atac.cds))
genome <- seqlengths(atac.blasts)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num = 100, window = 500000)
head(conns)
ccans <- generate_ccans(conns) 
head(ccans)

#saveRDS(ccans, file = "CCAN_Network.RDS")
#saveRDS(conns, file = "CCON_Scores.RDS")

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(atac.blasts) <- links
#saveRDS(atac.blasts, "ATAC_Blasts_1146_LINKS.RDS")

#Read in the links and do stuff with it 
atac.blasts <- readRDS("ATAC_Blasts_1146_LINKS.RDS")
Idents(atac.blasts) = "comparison.bmp.vs.t.specified"
table(atac.blasts$comparison.bmp.vs.t.specified)
diff_peaks <- FindMarkers(atac.blasts, ident.1 = "BMP-like-NR", ident.2 = "T-specified-R", 
                          only.pos = F)
write.table(diff_peaks, file = "Differential_Peaks.tsv", sep = '\t', quote = F)
diff_peaks$peak = rownames(diff_peaks)

#Analysis with all of BMP17 genes 
BMP17 = c("S100A4", "LGALS1", "FAM30A", "IGLL1", "MEF2C", 
          "CTSW", "PRDX1", "HCST", "HSH2D", "KLF2", "VAMP8", 
          "HOPX", "CYBA", "MT-ND3", "ENO1", "PEBP1", "CD44")
BMP_larger <- c("PRDX1", "S100A4", "BCL11A", "RAMP1", "SPINK2", "IGFBP7", "CD44", "C1QTNF4", "DNAJC1", "COMMD3", 
                              "PRSS57", "IGLL1", "LGALS1", "GYPC", "HCST", "ITGA4", "CTSW", "KIAA0087", "CYBA", "MEF2C", 
                              "MATK", "IFITM2", "CYTL1", "BAALC", "SMIM24", "HOXA10", "PTH2", "TENT5A", "CAPG", "PEBP1", 
                              "HSH2D", "VAMP8", "LPCAT2", "EEF1A2", "HOXA9", "FAM30A", "BMI1", "GOLGA8N", "DDAH2", "MT-ND3", 
                              "MT-ND4L", "MCTP2", "ENO1", "AHNAK", "TRH", "XIST", "CEBPG", "IGKC", "HOPX", "GZMA", "MSI2", 
                              "GSTP1", "MTRNR2L8", "APLP2", "PTPRE", "ZNF503", "AC103591.3", "KLF2", "DAD1", "CST3", "PLCG2", 
                              "TRPS1", "PRDX2", "ARL4C", "MDK", "ACY3")
T.spec_larger  <- c("HES4", "MAL", "SCGB3A1", "BCL11B", "BTG3", "TASP1", "CHI3L2", "TRGC2", "PMEPA1", "IFITM3",
                                  "NDFIP1", "SDCBP", "BEX3", "TCF7", "SH3BP5", "LAT", "MYO7B", "ITGA1", "CD3E", "HERPUD1",
                                  "CRYBG1", "PCGF5", "PMAIP1", "CD99", "MIR4435-2HG", "MZB1", "SIT1", "LY6H", "MGAT4A",
                                  "PRSS2", "L1TD1", "LCK", "ITM2B", "UBE2B", "CLDN5", "CCND3", "LPAR6", "SNTG2", "YBX3",
                                  "PCDH10", "LTB", "MALT1", "TSPYL2", "SAT1", "DNTT", "RBM38", "PPP1CB", "MAP1A", "SDF2L1",
                                  "DDIT4", "CYTOR", "TUBB2A", "HSP90B1")


geneconversions = AnnotationDbi::select(hs, keys = T.spec_larger, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
geneids = geneconversions$ENTREZID
mygenes.transcripts = subset(genes(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene), columns=c("tx_id", "tx_name","gene_id")), gene_id %in% geneids)
mygenes.tss = resize(mygenes.transcripts, width=1, fix='start')
expanded_bmp_granges <- resize(mygenes.tss, width=4000, fix="center")
print(as.data.frame(expanded_bmp_granges))

#Based on modules -----
links_table = as.data.frame(Links(atac.blasts))
ccan_all = readRDS("CCAN_Network.RDS")
split_coordinates <- strsplit(ccan_all$Peak, "-")
chromosome <- sapply(split_coordinates, "[[", 1)
start_pos <- as.integer(sapply(split_coordinates, "[[", 2))
end_pos <- as.integer(sapply(split_coordinates, "[[", 3))
ccan_gr <- GRanges(seqnames = chromosome, ranges = IRanges(start = start_pos, end = end_pos, 
                                                           CCAN = ccan_all$CCAN))
overlap_indices <- findOverlaps(ccan_gr, expanded_bmp_granges)

subset_coords <- ccan_gr[queryHits(overlap_indices), ]
print(as.data.frame(subset_coords))
all_groups <- unique(c(subset_coords$CCAN))

bmp_groups <- ccan_gr[ccan_gr$CCAN %in% all_groups, ] 
bmp_Granges <- makeGRangesFromDataFrame(bmp_groups)

peaks <- atac.blasts@assays$ATAC@ranges
peaks_df = as.data.frame(peaks)
bmp_Granges_df = as.data.frame(bmp_Granges)
overlaps <- findOverlaps(peaks, bmp_Granges)
main_peaks <- peaks[unique(overlaps@from)]
main_peaks_df = as.data.frame(main_peaks)
export.bed(main_peaks, con="TSpec_larger_CCAN.bed")

#Based on EP cutoff ------
ccon_scores = na.omit(readRDS("CCON_Scores.RDS"))
dim(ccon_scores)
sig_ccon = ccon_scores[ccon_scores$coaccess>=0.15, ]
dim(sig_ccon)

peaks <- atac.blasts@assays$ATAC@ranges
peaks_df = as.data.frame(peaks)
overlaps <- findOverlaps(peaks, expanded_bmp_granges)
peak_in_bmp <- peaks[queryHits(overlaps), ]
bmp_peak_df <- as.data.frame(peak_in_bmp)
bmp_peak_df$full = NULL
for(i in 1:nrow(bmp_peak_df)){
  bmp_peak_df$full[i] = paste(bmp_peak_df$seqnames[i], bmp_peak_df$start[i], bmp_peak_df$end[i], sep = "-")
}
ccon_bmp = sig_ccon[sig_ccon$Peak1 %in% bmp_peak_df$full, ]
dim(ccon_bmp)
bmp_ep_peaks = unique(c(as.vector(ccon_bmp$Peak1), as.vector(ccon_bmp$Peak2)))

split_coordinates <- strsplit(bmp_ep_peaks, "-")
chromosome <- sapply(split_coordinates, "[[", 1)
start_pos <- as.integer(sapply(split_coordinates, "[[", 2))
end_pos <- as.integer(sapply(split_coordinates, "[[", 3))
bmp_ep_peaks_gr <- GRanges(seqnames = chromosome, ranges = IRanges(start = start_pos, end = end_pos))
export.bed(bmp_ep_peaks_gr, con="TSpec_larger_EP_0.15.bed")

#Add fragments 
seurat.atac.frag = readRDS("/mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/Scripts/SeuratObj/Patients/Signac_30ETPsample_downsample_forPlot.rds")
table(seurat.atac.frag$comparison.bmp.vs.t.specified)
table(seurat.atac.frag$ctype_patient)

table(seurat.atac.frag$ctype_patient)
table(seurat.atac.frag$clinical.outcome)
table(seurat.atac.frag$leuk.subtype)
table(seurat.atac.frag$projected_ctype)

seurat.atac.frag@meta.data <- seurat.atac.frag@meta.data %>% 
  mutate(comparison.bmp.vs.t.specified = case_when(ctype_patient == "Blast" & D29.MRD > 0.01 & projected_ctype %in% c("HSPC", "LMPP", "CLP", "ETP") & leuk.subtype == "ETP T-ALL" ~ "BMP-like-NR",
                                                   ctype_patient == "Blast" & D29.MRD == 0 & projected_ctype %in% c("Pro-T", "Pre-T") & leuk.subtype == "ETP T-ALL" ~ "T-specified-R"))
seurat.atac.frag@meta.data$comparison.bmp.vs.t.specified[is.na(seurat.atac.frag$comparison.bmp.vs.t.specified)] = "Other"
table(seurat.atac.frag$comparison.bmp.vs.t.specified)

Links(seurat.atac.frag) <- Links(atac.blasts)

frags = Fragments(seurat.atac.frag)

frags = c()
samples = unique(seurat.atac.frag$sample)
for(item in samples){
  print(item)
  object = subset(seurat.atac.frag, subset = sample == as.character(item))
  frags[[item]] = CreateFragmentObject(path = paste0("ATAC_Fragments/T_ALL_",as.character(item),"_scATAC.fragments.downsample.txt.sort.bed.gz"), 
                               cells = colnames(object), validate.fragments = T)
}
Fragments(seurat.atac.frag) = NULL
Fragments(seurat.atac.frag) = frags

saveRDS(seurat.atac.frag, "Downsampled_with_Fragments.RDS")
seurat.atac.frag.blasts = subset(seurat.atac.frag, subset = ctype_patient == "Blast")
saveRDS(seurat.atac.frag.blasts, "Downsampled_Blasts_with_Fragments.RDS")
table(seurat.atac.frag.blasts$comparison.bmp.vs.t.specified)
Links(seurat.atac.frag.blasts) <- Links(atac.blasts)

#saveRDS(seurat.atac.frag.blasts, "Downsampled_Blasts_with_Fragments.RDS")

bmp_tspec = subset(seurat.atac.frag, subset = comparison.bmp.vs.t.specified == "Other", invert = T)
Links(bmp_tspec) = Links(seurat.atac.frag.blasts)

Idents(bmp_tspec) = 'comparison.bmp.vs.t.specified'

#Make figures 
diff_peaks = read.table("Differential_Peaks.tsv", sep = '\t', header = T)
diff_peaks$peaks = rownames(diff_peaks)

seurat.atac.frag.blasts = readRDS("Downsampled_Blasts_with_Fragments.RDS")
bmp_tspec = subset(seurat.atac.frag.blasts, subset = comparison.bmp.vs.t.specified == "Other", invert = T)
Links(bmp_tspec) = Links(seurat.atac.frag.blasts)
Idents(bmp_tspec) = 'comparison.bmp.vs.t.specified'
cov_plot = CoveragePlot(object=bmp_tspec, ymax = 500, 
                        region = "CD44", annotation = TRUE, group.by = 'comparison.bmp.vs.t.specified',
                        peaks = TRUE, title = TRUE, links = F, 
                        extend.upstream = 100000, extend.downstream = 0)
link_plot <- LinkPlot(
  object = bmp_tspec, min.cutoff = 0.15, 
  region = "CD44", extend.upstream = 100000, extend.downstream = 0
)
p1 = CombineTracks(list(cov_plot, link_plot))
p1
ggsave("ATAC_Figures/CD44Coverage_Plot_500.pdf", plot = p1, device = 'pdf', height = 4.5, width = 9)

data <- read.delim("BMP_larger_CCAN_output/knownResults.txt", header = TRUE, stringsAsFactors = FALSE)
data$FDR = p.adjust(exp(data$Log.P.value), method = "BH")
motifs <- sub("\\s*\\(.*", "", data$Motif.Name)
q_values <- data$FDR
motif_data <- data.frame(Motif = motifs, QValue = q_values)

motif_data <- motif_data[order(motif_data$QValue), ]
motif_data <- motif_data[motif_data$QValue < 0.01, ]
motif_data$QValue <- -log10(motif_data$QValue)
motif_data$Motif <- factor(
  motif_data$Motif, 
  levels = motif_data$Motif[order(motif_data$QValue)]
)
p1 <- ggplot(motif_data, aes(x = Motif, y = QValue)) +
  geom_bar(stat = "identity", fill = "darkred") +
  labs(x = "Motif", y = "-log10(FDR)") +
  theme_bw() + coord_flip()
p1
ggsave("ATAC_Figures/BMP_CCAN_Cicero_Motifs.pdf", plot = p1, device = "pdf", height = 4, width = 4)

