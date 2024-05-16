
setwd("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/8_remakeUMAP")
#make umap of bulkRNAseq
bulk.rna.object <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/1_harmonizemtx_make_seurat/vst.no.relapse.bulk.rna.umap.20.rds")
bulk.rna.object@reductions

bulk.rna.object = bulk.rna.object %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 1000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20, min.dist = 0.5)


bulk.rna.object = bulk.rna.object %>% 
  RunUMAP(dims = 1:10, min.dist = 0.5)
ElbowPlot(bulk.rna.object, ndims = 50)

  
DimPlot(bulk.rna.object, group.by = "ETP") + NoAxes() + NoLegend()
DimPlot(bulk.rna.object, group.by = "sc.data")
DimPlot(bulk.rna.object, group.by = "In.TARGET.Cohort..n..264..RNASeq.and.WES")
FeaturePlot(bulk.rna.object, "HOXA9")

coords.bulk.rna = bulk.rna.object@meta.data %>% cbind(bulk.rna.object@reductions$umap@cell.embeddings)


bulk.umaps.save = list()
bulk.umaps.save$shapes = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 1, aes(shape = ETP), color = ifelse(coords.bulk.rna$ETP == "ETP", yes = "#BE433D", no = "gray"))+
  geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T), size = 3.5, alpha = 1) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74", "Gray")) +
  theme(legend.position = "bottom")+
  scale_shape_manual(values = c(3, 2, 0, 0))




bulk.umaps.save$unknown = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 1, aes(shape = ETP), color = ifelse(coords.bulk.rna$ETP == "Unknown", yes = "#BE433D", no = "gray"))+
  #geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T), size = 4, alpha = 1) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74", "Gray")) +
  theme(legend.position = "bottom")+
  scale_shape_manual(values = c(3, 2, 0, 0))

bulk.umaps.save$subtype = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 1, aes(shape = ETP, color = ETP))+
  #geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T), size = 4, alpha = 1) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "black", "gray")) +
  theme(legend.position = "bottom")+
  scale_shape_manual(values = c(3, 2, 0, 0))

bulk.umaps.save$TARGET = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 1, aes(shape = ETP, color = In.TARGET.Cohort..n..264..RNASeq.and.WES))+
  #geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T), size = 4, alpha = 1) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("gray", "navy blue")) +
  theme(legend.position = "bottom")+
  scale_shape_manual(values = c(3, 2, 0, 0))

bulk.umaps.save$samplesource = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 1, aes(shape = ETP, color = Tumor.Specimen.Type))+
  #geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T), size = 4, alpha = 1) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("gray", "navy blue")) +
  theme(legend.position = "bottom")+
  scale_shape_manual(values = c(3, 2, 0, 0))


bulk.umaps.save$NearETP = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 1, aes(shape = ETP), color = ifelse(coords.bulk.rna$ETP == "Near-ETP", yes = "#DDA164", no = "gray"))+
  geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T), size = 3.5, alpha = 1) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "#B89B74", "gray")) +
  theme(legend.position = "bottom")+
  scale_shape_manual(values = c(3, 2, 0, 0))

bulk.umaps.save$NearETP2 = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 1, aes(shape = ETP,  color = ETP))+
  geom_point(aes(color = ETP), data = coords.bulk.rna %>% filter(sc.data == T), size = 3.5, alpha = 1) + 
  theme_classic() +
  #geom_density_2d(binwidth = 0.025)+
  scale_color_manual(values = c("#BE433D", "#DDA164", "gray", "gray")) +
  theme(legend.position = "bottom")+
  scale_shape_manual(values = c(3, 2, 0, 0))

getwd()
pdf("umap.bulk.rna2.pdf", height = 5, width = 4.5)
bulk.umaps.save
dev.off()

pdf("umap.bulk.rna.pdf_small.pdf", height = 5/2, width = 4.5/2)
bulk.umaps.save
dev.off()

t.all.oncogenes = c("TLX1", "TLX3", "NKX2-1", "TAL1", "TAL2", "LMO1", "LMO2", "LYL1", "FLT3", "HOXA9", "HOXA13")
t.all.onco.umap.bulk = FeaturePlot(bulk.rna.object, t.all.oncogenes, raster = F, order = T, combine = F)
colors = c("gray", (brewer.pal(n = 11, name = "Reds")))

for(i in 1:length(t.all.onco.umap.bulk)){
  t.all.onco.umap.bulk[[i]] = t.all.onco.umap.bulk[[i]] + NoAxes() + NoLegend() +
    scale_colour_gradientn(colours = colors)    
}

#library(ggrastr)
pdf("umap.bulk.rna.t.all.onco.pdf", width = 12, height =10)
cowplot::plot_grid(plotlist = t.all.onco.umap.bulk)
dev.off()




mutations = read_csv("/mnt/isilon/tan_lab/xuj5/PDX_Profiling/RScripts/06_select_PDX_by_ProgenitorProp/bm.prop.with.mutations.csv")

sc.samples = coords.bulk.rna %>% filter(sc.data == T) %>% left_join(mutations %>% dselect("USI" = sample.name, sample.group))
sc.samples$sample.name = sc.samples$USI
r.nr.umap = list()
r.nr.umap$all.umap = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna$ETP == "Non-ETP", yes = 0.7, no = 1.2), color = "gray") +
  scale_shape_manual(values = c(3, 2, 0, 0 ))+
  geom_point(aes(color = sample.group), data = sc.samples %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), size = 3, alpha = 0.9) + 
  theme_classic() +
  scale_color_manual(values = c("#BE433D","#459DA5")) +
  theme(legend.position = "none") +
  geom_text_repel(data = sc.samples %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

r.nr.umap$highlight.etp = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna$ETP == "ETP", yes = 1.2, no = 0.7), 
             color = ifelse(coords.bulk.rna$ETP == "ETP", yes = "#BE433D", no = "gray")) +
  scale_shape_manual(values = c(3, 2, 0, 0 ))+
  theme_classic() +
  theme(legend.position = "none")

r.nr.umap$highlight.near.etp = ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna$ETP == "Near-ETP", yes = 1.2, no = 0.7), 
             color = ifelse(coords.bulk.rna$ETP == "Near-ETP", yes = "#DDA164", no = "gray")) +
  scale_shape_manual(values = c(3, 2, 0, 0 ))+
  theme_classic() +
  theme(legend.position = "none")


r.nr.umap$subtype =  ggplot(coords.bulk.rna, aes(x = UMAP_1, y = UMAP_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  geom_point(alpha = 0.8, aes(shape = ETP), size = ifelse(coords.bulk.rna$ETP == "Non-ETP", yes = 0.7, no = 1.2), color = "gray") +
  scale_shape_manual(values = c(3, 2, 0, 0 ))+
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
coords.etp.bulk$sample.name = coords.etp.bulk$USI

sc.samples.etp = coords.etp.bulk %>% filter(sc.data == T) %>% left_join(mutations %>% dselect(sample.name, sample.group))


r.nr.umap$etp.umap10 = ggplot(coords.etp.bulk, aes(x = umap10_1, y = umap10_2)) + 
  geom_point(color = "gray") +
  #geom_density_2d(binwidth = 0.03)+
  #scale_shape_manual(values = c(3, 2, 0))+
  geom_point(aes(color = sample.group), 
             data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), 
             size = 3, alpha = 0.9) + 
  theme_classic() +
  scale_color_manual(values = c("#BE433D","#459DA5")) +
  theme(legend.position = "none") +
  geom_text_repel(data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

r.nr.umap$etp.umap10.marrow = ggplot(coords.etp.bulk, aes(x = umap10_1, y = umap10_2)) + 
  #geom_density_2d(binwidth = 0.03)+
  scale_shape_manual(values = c(1, 2, 4))+
  geom_point(aes(color = sample.group), 
             data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse", "ETP Low MRD no relapse")), size = 3, alpha = 0.9) + 
  geom_point(aes(shape = D29.BM)) +
  theme_classic() +
  scale_color_manual(values = c("#BE433D","#459DA5")) +
  #theme(legend.position = "none") +
  geom_text_repel(data = sc.samples.etp %>% filter(sample.group %in% c("ETP High MRD no relapse")), 
                  aes(label = sample.name), max.overlaps = 100, color = "black")

coords.etp.bulk = etp.bulk@meta.data %>% cbind(etp.bulk@reductions$umap15@cell.embeddings)
coords.etp.bulk$sample.name = coords.etp.bulk$USI

sc.samples.etp = coords.etp.bulk %>% filter(sc.data == T) %>% left_join(mutations %>% dselect(sample.name, sample.group))

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

