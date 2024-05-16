setwd("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/1_harmonizemtx_make_seurat")

load("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/TALL_X01_normalized_1362.Rdata")

#TPM

TPM = gexp.tpm

#VST - batch corrected, loosely filtered, vst transformed

VST = gexp.vst

select = dplyr::select

remove.ensembl.version = function(ensg.gene.names = NA){
  ensg.gene.names <- gsub('\\..+$', '', ensg.gene.names)
  return(ensg.gene.names)
}

move.to.rowname = function(df = gene.level.etp, column = "hgnc_symbol"){
  df = df %>% as.data.frame()
  rownames(df) = df[[column]]
  df = df %>% dplyr::select(-.data[[column]])
  return(df)
}



convert.petri.to.10x = function(raw.matrix = NA){ #input raw matrix, genes by samples
  
  
  print("reading in 10X GTF")
  #read in gtf
  single.cell.gtf = "/mnt/isilon/tan_lab/chenc6/Tools/SingleCellAnnotation/GRCh38/genes/genes.gtf"
  import.gtf = rtracklayer::import(single.cell.gtf)
  gtf.df = as.data.frame(import.gtf)
  
  #only genes - 33k
  print("filtering for genes in GTF")
  single.cell.genes.only.genes = gtf.df %>% filter(type == "gene")
  
  single.cell.genes.only.genes %>% dim()
  
  single.cell.genes.only.genes.conversion = single.cell.genes.only.genes %>% dplyr::select(gene_id, gene_name)
  
  
  #requires single cell.genes.only.genes.converstion, which is from 10x GTF
  
  #read in TPM / FPKM
  #raw.matrix = Txi_gene[["abundance"]] %>% as_tibble(rownames = "ENSG_Name")
  raw.matrix = raw.matrix
  print("loading raw matrix")
  print(raw.matrix[1:5,1:5])
  #rename column names to reflect samples
  
  
  
  
  #now split the gene names to reflect ensembl_id and st.jude gene name
  if(str_detect(rownames(raw.matrix)[1], pattern = "_")){
    print('naming first column of matrix to ENSG_Name')
    raw.matrix = raw.matrix %>% as_tibble(rownames = "ENSG_Name")
    raw.matrix = raw.matrix %>% mutate(ensembl.id = str_split_fixed(ENSG_Name, pattern = "_", n = 2)[,1],
                                       st.jude.gene.name = str_split_fixed(ENSG_Name, pattern = "_", n = 2)[,2])
    dictionary = raw.matrix %>% dplyr::select(ENSG_Name, ensembl.id, st.jude.gene.name)
    
  }else{
    
    raw.matrix = raw.matrix %>% as_tibble(rownames = "ensembl.id")
    dictionary = raw.matrix %>% dplyr::select(ensembl.id)
    
  }
  
  
  #make dictionary to convert with 10x genes
  dictionary = dictionary %>% mutate(ensembl.id.no.version = remove.ensembl.version(ensembl.id))
  dictionary$in.10X.ref = dictionary$ensembl.id.no.version %in% single.cell.genes.only.genes.conversion$gene_id
  table(dictionary$in.10X.ref) %>% print()#check how many genes are preserved 
  
  #add 10x gene annos based on ensembl.id.no.version
  dictionary = dictionary %>% left_join(single.cell.genes.only.genes.conversion %>% dplyr::select("ensembl.id.no.version" = gene_id, "Gene.Name.10X" = gene_name))
  table(is.na(dictionary$Gene.Name.10X))
  
  #filter dictionary to only have genes to keep
  dictionary.genes.to.keep = dictionary %>% filter(in.10X.ref == T)
  dictionary.genes.to.keep
  
  #filter harmonized matrix to only have our genes
  harmonized.matrix.tpm = raw.matrix[raw.matrix$ensembl.id %in% dictionary.genes.to.keep$ensembl.id,] #filter to our genes
  harmonized.matrix.tpm %>% dim()# now 33k by 1046
  harmonized.matrix.tpm = harmonized.matrix.tpm %>% left_join(dictionary.genes.to.keep %>% dplyr::select(ensembl.id, Gene.Name.10X)) #change to our gene names
  
  #move 10x gene name to front and only keep columns with values
  harmonized.matrix.tpm = harmonized.matrix.tpm %>% select(Gene.Name.10X, everything())
  harmonized.matrix.tpm = harmonized.matrix.tpm %>% select(-ensembl.id)
  
  #move.to.rowname(df =harmonized.matrix.tpm, column = "Gene.Name.10X") #doesnt work due to duplicated rownames
  
  #sum up duplicated rowanmes
  harmonized.matrix.tpm.duplicates.removed = dplyr::group_by(harmonized.matrix.tpm, Gene.Name.10X) %>% dplyr::summarise_all(sum)
  harmonized.matrix.tpm.duplicates.removed.mtx = harmonized.matrix.tpm.duplicates.removed %>% move.to.rowname("Gene.Name.10X") #move gene names to front
  
  #return matrix
  return(harmonized.matrix.tpm.duplicates.removed.mtx)
}

harmonized.matrix.tpm.duplicates.removed.mtx = convert.petri.to.10x(TPM)


harmonized.matrix.tpm.duplicates.removed.mtx[1:5, 1:5]

harmonized.matrix.vst.duplicates.removed.mtx = convert.petri.to.10x(VST)
harmonized.matrix.vst.duplicates.removed.mtx[1:5, 1:5]


#dir.save =getwd()
dir.save =  "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/1_harmonizemtx_make_seurat"
filename = paste0(dir.save, "/", nrow(harmonized.matrix.tpm.duplicates.removed.mtx), "_genes_", ncol(harmonized.matrix.tpm.duplicates.removed.mtx), "_samples_10XHarmonized_TPM_noLog.rds") %>% print()
saveRDS(harmonized.matrix.tpm.duplicates.removed.mtx, filename)


filename2 = paste0(dir.save, "/", nrow(harmonized.matrix.vst.duplicates.removed.mtx), "_genes_", ncol(harmonized.matrix.vst.duplicates.removed.mtx), "_samples_10XHarmonized_VST_noLog.rds") %>% print()
saveRDS(harmonized.matrix.vst.duplicates.removed.mtx, filename2)


colSums(harmonized.matrix.tpm.duplicates.removed.mtx) %>% head()
colSums(harmonized.matrix.vst.duplicates.removed.mtx) %>% head()
colSums(VST) %>% head()
colSums(TPM) %>% head()


#make seurat object
load("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/data_sharing/TALL_X01_patient_annotations_1362.Rdata")

patient_anno


final.cohort.clinical.genomics <- read_excel("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/AALL0434 X01 and ALSF Abbreviated Dataset December 2022 no PHI.xlsx", sheet = "Final Cohort w Genomics no PHI")
final.cohort.clinical.genomics$ETP.STATUS %>% table()



table(colnames(TPM) %in% eset_counts@phenoData@data$sample.name)
colnames(TPM)[!colnames(TPM) %in% eset_counts@phenoData@data$sample.name]
table(colnames(TPM) %in% patient_anno$USI)

table(colnames(TPM) %in% final.cohort.clinical.genomics$USI)
colnames(TPM)[!colnames(TPM) %in% final.cohort.clinical.genomics$USI]
colnames(TPM)[!colnames(TPM) %in% final.cohort.clinical.genomics$USI]


#go with patient_anno

clinical.t.all.40 = read_csv(files$clinical.mdata)

colnames(TPM) == patient_anno$USI

patient_anno_df = patient_anno %>% move.to.rowname("USI")
patient_anno_df$USI = patient_anno$USI
patient_anno_df$In.TARGET.Cohort..n..264..RNASeq.and.WES[patient_anno_df$In.TARGET.Cohort..n..264..RNASeq.and.WES %>% is.na()] = "No"
patient_anno_df$sc.data = patient_anno_df$USI %in% clinical.t.all.40$sample.name


bulk.rna.object = CreateSeuratObject(counts = harmonized.matrix.tpm.duplicates.removed.mtx, 
                                     assay = "bulkRNA_TPM", 
                                     meta.data = patient_anno_df)

bulk.rna.object@assays$bulkRNA_VST = CreateAssayObject(counts = harmonized.matrix.vst.duplicates.removed.mtx)


bulk.rna.object.TPM = CreateSeuratObject(counts = harmonized.matrix.tpm.duplicates.removed.mtx, 
                                     assay = "bulkRNA_TPM", 
                                     meta.data = patient_anno_df)

bulk.rna.object.TPM = bulk.rna.object.TPM %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 500) %>% 
  ScaleData(vars.to.regress = "In.TARGET.Cohort..n..264..RNASeq.and.WES") %>% 
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

table(bulk.rna.object.TPM$In.TARGET.Cohort..n..264..RNASeq.and.WES)

pdf("TPM_dataset_QC.pdf", height = 4, width = 12)
DimPlot(bulk.rna.object.TPM, group.by = "ETP.STATUS") +
DimPlot(bulk.rna.object.TPM, group.by = "In.TARGET.Cohort..n..264..RNASeq.and.WES")+
DimPlot(bulk.rna.object.TPM, group.by = "sc.data")
dev.off()

#FeaturePlot(bulk.rna.object, "HOXA9")


bulk.rna.object.vst = CreateSeuratObject(counts = harmonized.matrix.vst.duplicates.removed.mtx, 
                                     assay = "bulkRNA_VST", 
                                     meta.data = patient_anno_df)

bulk.rna.object.vst = bulk.rna.object.vst %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 500) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

pdf("VST_dataset_QC.pdf",height = 4, width = 12)
DimPlot(bulk.rna.object.vst, group.by = "ETP.STATUS") +
DimPlot(bulk.rna.object.vst, group.by = "In.TARGET.Cohort..n..264..RNASeq.and.WES")+
DimPlot(bulk.rna.object.vst, group.by = "sc.data")
dev.off()


#save RDS files of TPM and VST matrices
saveRDS(bulk.rna.object.vst, "vst.bulk.rna.umap.20.rds")
saveRDS(bulk.rna.object.TPM, "TPM.bulk.rna.umap.20.rds")
saveRDS(bulk.rna.object, "TPM.vst.bulk.rna.umap.20.rds")


saveRDS(bulk.rna.object.vst, "vst.bulk.rna.umap.20.rds")

#exclude relapse
bulk.rna.object.vst$USI %>% grep(pattern = "_R")
bulk.rna.object.vst$relapse.primary = "primary"
bulk.rna.object.vst$relapse.primary[grep(x = bulk.rna.object.vst$USI, pattern = "_R")] = "relapse"
table(bulk.rna.object.vst$relapse.primary)




bulk.rna.object.vst.no.relapse = subset(bulk.rna.object.vst, relapse.primary == "primary")
bulk.rna.object.vst.no.relapse = bulk.rna.object.vst.no.relapse %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 500) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

colnames(bulk.rna.object.vst.no.relapse)[!colnames(bulk.rna.object.vst.no.relapse) %in% final.cohort.clinical.genomics$USI]
bulk.rna.object.vst.no.relapse@meta.data %>% filter(USI %in% colnames(bulk.rna.object.vst.no.relapse)[!colnames(bulk.rna.object.vst.no.relapse) %in% final.cohort.clinical.genomics$USI])

etp.status = bulk.rna.object.vst.no.relapse@meta.data %>% dselect(USI, ETP.STATUS) %>% 
  left_join(final.cohort.clinical.genomics %>% dselect(USI, "teachey.ETP.STATUS" = ETP.STATUS, ETP_status)) %>% 
  left_join(eset_counts@phenoData@data %>% dselect("USI" = sample.name, "old.ETP" = ETP))

etp.status%>% filter(USI %in% colnames(bulk.rna.object.vst.no.relapse)[!colnames(bulk.rna.object.vst.no.relapse) %in% final.cohort.clinical.genomics$USI])

table(etp.status$ETP.STATUS == etp.status$teachey.ETP.STATUS)

pdf("Norelapse_VST_dataset_QC.pdf",height = 4, width = 12)
DimPlot(bulk.rna.object.vst.no.relapse, group.by = "ETP.STATUS") +
  DimPlot(bulk.rna.object.vst.no.relapse, group.by = "In.TARGET.Cohort..n..264..RNASeq.and.WES")+
  DimPlot(bulk.rna.object.vst.no.relapse, group.by = "sc.data")
dev.off()

saveRDS(bulk.rna.object.vst.no.relapse, "vst.no.relapse.bulk.rna.umap.20.rds")

etp.bulk.vst = subset(bulk.rna.object.vst.no.relapse, ETP.STATUS == "ETP")

saveRDS(etp.bulk.vst, "etp.110.bulk.vst.rds")










#make umap of bulkRNAseq
eset_counts <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/bulkRNA/13_FinalMtx/eset_counts.rds")
counts_15732_genes_1373_samples_10XHarmonized_TPM_noLog <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/bulkRNA/13_FinalMtx/15732_genes_1373_samples_10XHarmonized_TPM_noLog.rds")



bulk.rna.object = CreateSeuratObject(counts = counts_15732_genes_1373_samples_10XHarmonized_TPM_noLog, 
                                     assay = "bulkRNA", 
                                     meta.data = eset_counts@phenoData@data)

bulk.rna.object = bulk.rna.object %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 500) %>% 
  ScaleData(vars.to.regress = "in.target") %>% 
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

DimPlot(bulk.rna.object, group.by = "ETP")
DimPlot(bulk.rna.object, group.by = "sc.data")
DimPlot(bulk.rna.object, group.by = "in.target")
FeaturePlot(bulk.rna.object, "HOXA9")

coords.bulk.rna = bulk.rna.object@meta.data %>% cbind(bulk.rna.object@reductions$umap@cell.embeddings)
umap.fig.1b = list()


umap.fig.1b$plain = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.3)+
  geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.8) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74"))

umap.fig.1b$nolabel = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.3)+
  geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.8) + 
  theme_classic() +
  geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74"))



umap.fig.1b$labeled.etp.mrd.1= ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2, color = ETP)) + 
  geom_point(alpha = 0.3)+
  geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3) + 
  geom_density_2d(binwidth = 0.03)+
  theme_classic() +
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74")) +
  geom_text_repel(data = coords.bulk.rna %>% filter(sc.data == T, relapse.rna.seq == "Primary", ETP == "ETP", D29.MRD > 1), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")



fig1.save = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure1/"
pdf(paste0(fig1.save, "umap.bulk.rna.pdf"))
umap.fig.1b
dev.off()


umap.fig.1b$shapes = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.5, shape = 4)+
  geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.8) + 
  theme_classic() +
  geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74")) +
  theme(legend.position = "none")

pdf(paste0(fig1.save, "umap.bulk.rna.pdf"))
umap.fig.1b
dev.off()



bulk.rna.object = bulk.rna.object %>% RunUMAP(dims = 1:20, min.dist = 0.5, reduction.name = "umap20_mindist05")
DimPlot(bulk.rna.object, group.by = "ETP", reduction = "umap20_mindist05")

coords.bulk.rna.min.dist.05 = bulk.rna.object@meta.data %>% cbind(bulk.rna.object@reductions$umap20_mindist05@cell.embeddings)
umap.fig.1b.min.dist.05 = list()


umap.fig.1b.min.dist.05$plain = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.3)+
  geom_point(aes(color = ETP), data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.8) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74"))

umap.fig.1b.min.dist.05$nolabel = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.3)+
  geom_point(aes(color = ETP), data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.8) + 
  theme_classic() +
  geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74"))



umap.fig.1b.min.dist.05$labeled.etp.mrd.1= ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2, color = ETP)) + 
  geom_point(alpha = 0.3)+
  geom_point(aes(color = ETP), data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3) + 
  geom_density_2d(binwidth = 0.03)+
  theme_classic() +
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74")) +
  geom_text_repel(data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary", ETP == "ETP", D29.MRD > 1), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")


umap.fig.1b.min.dist.05$shapes = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.5, shape = 4)+
  geom_point(aes(color = ETP), data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.8) + 
  theme_classic() +
  geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74")) +
  theme(legend.position = "none")


ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna.min.dist.05$ETP == "Non-ETP", yes = 0.7, no = 1.2)) +
  scale_shape_manual(values = c(3, 9, 0))+
  geom_point(aes(color = ETP), data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.9) + 
  theme_classic() +
  geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74")) +
  theme(legend.position = "none")


pdf(paste0(fig1.save, "umap.bulk.rna.mindist05.pdf"))
umap.fig.1b.min.dist.05
dev.off()

final.1b = list()
final.1b$color = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna.min.dist.05$ETP == "Non-ETP", yes = 0.7, no = 1.2)) +
  scale_shape_manual(values = c(3, 9, 0))+
  geom_point(aes(color = ETP), data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.9) + 
  theme_classic() +
  geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74")) +
  theme(legend.position = "none")

final.1b$gray = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna.min.dist.05$ETP == "Non-ETP", yes = 0.7, no = 1.2), color = "gray") +
  #scale_shape_manual(values = c(3, 9, 0))+
  scale_shape_manual(values = c(3, 2, 0))+
  geom_point(aes(color = ETP), data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.9) + 
  theme_classic() +
  geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74")) +
  theme(legend.position = "none")

final.1b$label = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2, color = ETP)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna.min.dist.05$ETP == "Non-ETP", yes = 0.7, no = 1.2), color = "gray") +
  scale_shape_manual(values = c(3, 2, 0))+
  geom_point(aes(color = ETP), data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary"), size = 3, alpha = 0.9) + 
  theme_classic() +
  geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74")) +
  theme(legend.position = "none") +
  geom_text_repel(data = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T, relapse.rna.seq == "Primary", ETP == "ETP", D29.MRD > 1), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

pdf(paste0(fig1.save, "shape_umap.bulk.rna.mindist05.pdf"), height = 4, width = 4)
final.1b
dev.off()

t.all.oncogenes = c("TLX1", "TLX3", "NKX2-1", "TAL1", "TAL2", "LMO1", "LMO2", "LYL1", "FLT3", "HOXA9", "HOXA13")
t.all.onco.umap.bulk = FeaturePlot(bulk.rna.object, t.all.oncogenes, raster = T)

pdf(paste0(fig1.save, "umap.bulk.rna.t.all.onco.pdf"), width = 10, height =10)
t.all.onco.umap.bulk
dev.off()

bulk.rna.object = bulk.rna.object %>% ScaleData(features = t.all.oncogenes)
t.all.onco.umap.bulk.mindist0.5 = FeaturePlot(bulk.rna.object, t.all.oncogenes, raster = T, reduction = "umap20_mindist05", slot = "scale.data")
for(i in 1:length(t.all.onco.umap.bulk.mindist0.5)){
  t.all.onco.umap.bulk.mindist0.5[[i]] = t.all.onco.umap.bulk.mindist0.5[[i]] + NoAxes()
}

pdf(paste0(fig1.save, "umap.bulk.rna.t.all.onco.mindist05.pdf"), width = 10, height =10)
t.all.onco.umap.bulk.mindist0.5
dev.off()


mutations = read_csv("/mnt/isilon/tan_lab/xuj5/PDX_Profiling/RScripts/06_select_PDX_by_ProgenitorProp/bm.prop.with.mutations.csv")

sc.samples = coords.bulk.rna.min.dist.05 %>% filter(sc.data == T,relapse.rna.seq == "Primary") %>% left_join(mutations)

r.nr.umap = list()
r.nr.umap$all.umap = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna.min.dist.05$ETP == "Non-ETP", yes = 0.7, no = 1.2), color = "gray") +
  scale_shape_manual(values = c(3, 2, 0))+
  geom_point(aes(color = sample.group), data = sc.samples %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), size = 3, alpha = 0.9) + 
  theme_classic() +
  scale_color_manual(values = c("#BE433D","#459DA5")) +
  theme(legend.position = "none") +
  geom_text_repel(data = sc.samples %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

r.nr.umap$highlight.etp = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna.min.dist.05$ETP == "ETP", yes = 1.2, no = 0.7), 
             color = ifelse(coords.bulk.rna.min.dist.05$ETP == "ETP", yes = "#BE433D", no = "gray")) +
  scale_shape_manual(values = c(3, 2, 0))+
  theme_classic() +
  theme(legend.position = "none")

r.nr.umap$highlight.near.etp = ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna.min.dist.05$ETP == "Near-ETP", yes = 1.2, no = 0.7), 
             color = ifelse(coords.bulk.rna.min.dist.05$ETP == "Near-ETP", yes = "#DDA164", no = "gray")) +
  scale_shape_manual(values = c(3, 2, 0))+
  theme_classic() +
  theme(legend.position = "none")


r.nr.umap$subtype =  ggplot(coords.bulk.rna.min.dist.05, aes(x = umap20_mindist05_1, y = umap20_mindist05_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna.min.dist.05$ETP == "Non-ETP", yes = 0.7, no = 1.2), color = "gray") +
  scale_shape_manual(values = c(3, 2, 0))+
  geom_point(aes(color = subtype), data = sc.samples %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), size = 3, alpha = 0.9) + 
  theme_classic() +
  #theme(legend.position = "none") +
  geom_text_repel(data = sc.samples %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")


etp.bulk = subset(bulk.rna.object, ETP == "ETP")


etp.bulk = etp.bulk%>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 1000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

etp.bulk = etp.bulk %>% RunUMAP(dims = 1:10, reduction.name = "umap10")
etp.bulk = etp.bulk %>% RunUMAP(dims = 1:15, reduction.name = "umap15")

coords.etp.bulk = etp.bulk@meta.data %>% cbind(etp.bulk@reductions$umap10@cell.embeddings)
sc.samples.etp = coords.etp.bulk %>% filter(sc.data == T,relapse.rna.seq == "Primary") %>% left_join(mutations)


r.nr.umap$etp.umap10 = ggplot(coords.etp.bulk, aes(x = umap10_1, y = umap10_2)) + 
  geom_point(color = "gray") +
  #geom_density_2d(binwidth = 0.03)+
  #scale_shape_manual(values = c(3, 2, 0))+
  geom_point(aes(color = sample.group), data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), size = 3, alpha = 0.9) + 
  theme_classic() +
  scale_color_manual(values = c("#BE433D","#459DA5")) +
  theme(legend.position = "none") +
  geom_text_repel(data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

r.nr.umap$etp.umap10.marrow = ggplot(coords.etp.bulk, aes(x = umap10_1, y = umap10_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  scale_shape_manual(values = c(1, 2, 4))+
  geom_point(aes(color = sample.group), data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), size = 3, alpha = 0.9) + 
  geom_point(aes(shape = D29.BM)) +
  theme_classic() +
  scale_color_manual(values = c("#BE433D","#459DA5")) +
  #theme(legend.position = "none") +
  geom_text_repel(data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

coords.etp.bulk = etp.bulk@meta.data %>% cbind(etp.bulk@reductions$umap15@cell.embeddings)
sc.samples.etp = coords.etp.bulk %>% filter(sc.data == T,relapse.rna.seq == "Primary") %>% left_join(mutations)

r.nr.umap$etp.umap15 =  ggplot(coords.etp.bulk, aes(x = umap15_1, y = umap15_2)) + 
  geom_point(color = "gray") +
  #geom_density_2d(binwidth = 0.03)+
  #scale_shape_manual(values = c(3, 2, 0))+
  geom_point(aes(color = sample.group), data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), size = 3, alpha = 0.9) + 
  theme_classic() +
  scale_color_manual(values = c("#BE433D","#459DA5")) +
  theme(legend.position = "none") +
  geom_text_repel(data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

r.nr.umap$etp.umap15.marrow = ggplot(coords.etp.bulk, aes(x = umap15_1, y = umap15_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  scale_shape_manual(values = c(1, 2, 4))+
  geom_point(aes(color = sample.group), data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), size = 3, alpha = 0.9) + 
  geom_point(aes(shape = D29.BM)) +
  theme_classic() +
  scale_color_manual(values = c("#BE433D","#459DA5")) +
  #theme(legend.position = "none") +
  geom_text_repel(data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

etp.bulk = etp.bulk %>% ScaleData(features = t.all.oncogenes)
etp.onco.umap.bulk.mindist0.5 = FeaturePlot(etp.bulk, t.all.oncogenes, raster = T, reduction = "umap", slot = "scale.data")
for(i in 1:length(t.all.onco.umap.bulk.mindist0.5)){
  t.all.onco.umap.bulk.mindist0.5[[i]] = t.all.onco.umap.bulk.mindist0.5[[i]] + NoAxes()
}


pdf(paste0("umap.bulk.rna.etp.pdf"), width = 5, height =5)
r.nr.umap
r.nr.umap$etp.umap15.marrow + theme(legend.position = "none")
r.nr.umap$etp.umap10.marrow + theme(legend.position = "none")
dev.off()


etp.bulk$D29.BM %>% table()






#save objects
setwd("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure1")
saveRDS(bulk.rna.object, "bulk.rna.object.umap20.rds")
saveRDS(etp.bulk, "etp.bulk.umap10.15.rds")