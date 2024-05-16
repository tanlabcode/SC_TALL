source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/JX_scFunctions.R')
setwd("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/5_BMP17_VST_AUC")

bulk.rna.object.umap20 <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/1_harmonizemtx_make_seurat/vst.no.relapse.bulk.rna.umap.20.rds")



expr.matrix = bulk.rna.object.umap20@assays$bulkRNA_VST@counts


#rank cells
cells_rankings_vst <- AUCell_buildRankings(expr.matrix, nCores = 12)
getwd()
saveRDS(cells_rankings_vst, file="cells_rankings_vst.rds")


#construct gene lists
geneSets = list()

#BMP17
ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9 <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9.rds")
consensus.t.spec.sig.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/consensus.t.spec.sig.nofilter.rds")
t.spec.10 = consensus.t.spec.sig.nofilter %>% filter(avg.fc > 0.9)

#BMP surface
consensus.bmp.surface.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/BMP_surface_phenotype/consensus.bmp.surface.nofilter.rds")

#Tspec Surface
consensus.t.spec.surface.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/BMP_surface_phenotype/consensus.t.spec.surface.nofilter.rds")

#BMP119
BMP.119.ETP <- read_csv("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/Figure2/119_gene_deg.sig.csv")

#TCF7LEF1

t.all.40.eps=read_tsv("/mnt/isilon/tan_lab/yuw1/R_work_dir/ETP_ALL/EP_Prediction2/regrRes4_EP_overall.txt")
edges = read_tsv("/mnt/isilon/tan_lab/yuw1/R_work_dir/ETP_ALL/EP_Prediction2/TRN_edges_enrichedTFs_in_Near-ETP_stringentPlus.txt")
nodes = read_tsv("/mnt/isilon/tan_lab/yuw1/R_work_dir/ETP_ALL/EP_Prediction2/TRN_nodes_enrichedTFs_in_Near-ETP_stringentPlus.txt")
ep.predictions.t.all.40 <- readRDS("/mnt/isilon/tan_lab/yuw1/R_work_dir/ETP_ALL/EP_Prediction2/regrRes4ep_prediction.rds")

t.all.40.eps %>% filter(gene_name == "CD5") %>% arrange(fdr)

edges %>% filter(gene_name == "CD5")


tcf7.lef1.targets = edges %>% filter(TF %in% c("TCF7", "LEF1"))
tcf7.lef1.targets.up = tcf7.lef1.targets %>% filter(direction == "up")
tcf7.lef1.targets.down = tcf7.lef1.targets %>% filter(direction == "down")

tcf7.lef1.cd5.activtors = edges %>% filter(gene_name %in% c("TCF7", "LEF1", "CD5")) 
tcf7.lef1.activtors = edges %>% filter(gene_name %in% c("TCF7", "LEF1")) 

c(tcf7.lef1.targets.up$gene_name, tcf7.lef1.cd5.activtors$TF) %>% unique()
tcf7.lef1.targets.down$gene_name %>% unique()


##FINAL GENE LISTS
geneSets$BMP17 = ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9$gene
geneSets$tspec10 = t.spec.10$gene
geneSets$BMPSurface.pos = consensus.bmp.surface.nofilter$coding
geneSets$BMPSurface.neg = consensus.t.spec.surface.nofilter$coding
geneSets$BMP.119.ETP.pos = BMP.119.ETP$gene[BMP.119.ETP$enriched == "BMP-like-NR"]
geneSets$BMP.119.ETP.neg = BMP.119.ETP$gene[BMP.119.ETP$enriched == "T-specified-R"]
geneSets$TCF7LEF1.pos =  c(tcf7.lef1.targets.up$gene_name, tcf7.lef1.cd5.activtors$TF) %>% unique()
geneSets$TCF7LEF1.neg =  tcf7.lef1.targets.down$gene_name %>% unique()

#Calculate AUC
cells_AUC.BMP <- AUCell_calcAUC(geneSets, cells_rankings_vst, 
                                aucMaxRank=nrow(cells_rankings_vst)*0.25, verbose = T)


bulk.rna.object.umap20[["AUCell_BMP"]] <- CreateAssayObject(cells_AUC.BMP@assays@data$AUC)


bulk.rna.object.umap20[["AUCell_BMP"]]@counts["BMPSurface.pos",] %>% summary()
bulk.rna.object.umap20[["AUCell_BMP"]]@counts["BMPSurface.neg",] %>% summary()

bulk.rna.object.umap20$BMP.surface.AUC = bulk.rna.object.umap20[["AUCell_BMP"]]@counts["BMPSurface.pos",] - bulk.rna.object.umap20[["AUCell_BMP"]]@counts["BMPSurface.neg",]
bulk.rna.object.umap20$BMP.119.ETP.AUC = bulk.rna.object.umap20[["AUCell_BMP"]]@counts["BMP.119.ETP.pos",] - bulk.rna.object.umap20[["AUCell_BMP"]]@counts["BMP.119.ETP.neg",]
bulk.rna.object.umap20$TCF7LEF1.AUC = bulk.rna.object.umap20[["AUCell_BMP"]]@counts["TCF7LEF1.pos",] - bulk.rna.object.umap20[["AUCell_BMP"]]@counts["TCF7LEF1.neg",]



#make list of bulk seurats
bulk.rna.object.umap20@meta.data = bulk.rna.object.umap20@meta.data %>% survival.conversion.for.plotting()
saveRDS(bulk.rna.object.umap20, "bulk.rna.object.umap20.with.BMP.scoring.rds")

metadata = cbind(bulk.rna.object.umap20@meta.data, bulk.rna.object.umap20[["AUCell_BMP"]]@counts %>% t())
metadata$sample.name = rownames(metadata)
metadata.for.haley = metadata %>% select(sample.name, BMP17, BMP.surface.AUC, tspec10, ETP)

bmp.high.cases = c()
for(status in c("ETP", "Near-ETP", "Non-ETP", "Unknown")){
  print(status)
  to.loop = metadata.for.haley %>% filter(ETP == status)
  BMP.high.cutoff =  quantile(to.loop$BMP17, probs = 0.70)
  bmp.high.cases = c(bmp.high.cases, to.loop$sample.name[to.loop$BMP17 >= BMP.high.cutoff])

}

metadata.for.haley$BMP.like.high = ifelse(metadata.for.haley$sample.name %in% bmp.high.cases, yes = "BMP.High", no = "BMP.low")
table(metadata.for.haley$BMP.like.high )

write_csv(metadata.for.haley, "primary-rna-seq-cases-with-BMP-sig.csv")

ggplot(metadata, aes(x =ETP, y = BMP17))+geom_violin()

#0434 ETP
etp.bulk = subset(bulk.rna.object.umap20, ETP == "ETP")

#0434 Near-ETP
near.etp.bulk = subset(bulk.rna.object.umap20, ETP == "Near-ETP")

#0434 Non-ETP
non.etp.bulk = subset(bulk.rna.object.umap20, ETP == "Non-ETP")

#0434 Unknown
unknown.bulk = subset(bulk.rna.object.umap20, ETP == "Unknown")


#0434 Non-ETP and Unknown
non.etp.plus.unknown.bulk = subset(bulk.rna.object.umap20, ETP %in% c("Non-ETP", "Unknown"))

bulk.rna.objects = list(bulk.rna.object.umap20, etp.bulk, near.etp.bulk, non.etp.bulk, unknown.bulk, non.etp.plus.unknown.bulk)
names(bulk.rna.objects) = c("0434 All", "0434 ETP", "0434 Near-ETP", "0434 Non-ETP", "0434 Unknown", "0434 Non-ETP-Unknown")


#do survival analysis


do.survival.given.gene.set = function(score.name, seurat.object, seurat.obj.name, lower.limit = 0.6){
  
  print('fixing surv metadata')
  
  #seurat.object@meta.data = seurat.object@meta.data %>% survival.conversion.for.plotting()
  seurat.object@meta.data = seurat.object@meta.data %>% mutate(status.OS.plot = case_when(status.OS == 0 ~ "Alive",
                                                                                          status.OS == 1 ~ "Dead"))
  
  
  
  #survival plot
  library(survival)
  library(ggfortify)
  
  
  
  print('running cox')
  survival.data = seurat.object@meta.data
  
  if(score.name %in% rownames(seurat.object@assays$AUCell_BMP@counts)){
    survival.data[[score.name]] = seurat.object@assays$AUCell_BMP@counts[score.name,]
  }
  
  
  survival.data$EOC.MRD = survival.data$EOC.MRD %>% as.numeric()
  
  survival.data$mpp.rank =  survival.data[[score.name]] %>% rank()
  survival.data$status.OS.plot %>% table()
  
  
  survival.data = survival.data %>% mutate(high.mpp = case_when(mpp.rank > (nrow(survival.data)*2/3) ~ "BMP High",
                                                                mpp.rank <= (nrow(survival.data)*1/3) ~ "BMP Low"))
  
  
  table(survival.data$high.mpp)
  survival.data$status.OS = as.numeric(survival.data$status.OS)
  os_km_fit <- survfit(Surv(time.OS, status.OS) ~ high.mpp , data=survival.data)
  cox.test = coxph(Surv(time.OS, status.OS) ~ high.mpp + D29.MRD, data=survival.data)
  cox.coefs.gene.and.mrd = summary(cox.test)
  
  if(nrow(survival.data) >100 ){
    p.val = cox.coefs.gene.and.mrd$logtest[3] %>% print()
  }else{
    p.val = cox.coefs.gene.and.mrd$sctest[3] %>% print()
    
  }
  
  
  
  if("cns.status" %in% colnames(survival.data)){
    cox.test2 = coxph(Surv(time.OS, status.OS) ~ survival.data[[score.name]] +  D29.MRD + cns.status+ dx.age + dx.wbc, data=survival.data)
    cox.test.summary = summary(cox.test2)
    p.val.signature = cox.test.summary$coefficients[1,5] %>% print()
    
  }else{
    cox.test2 = coxph(Surv(time.OS, status.OS) ~ survival.data[[score.name]] + D29.MRD, data=survival.data)
    cox.test.summary = summary(cox.test2)
    p.val.signature = cox.test.summary$coefficients[1,5] %>% print()
    
  }
  
  
  #ccp.val = cox.coefs.gene.and.mrd$sctest[3] %>% print()
  
  p.os = autoplot(os_km_fit)+ 
    labs(x = "OS (Days) ", y = "Survival Probability", 
         title = "") + 
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_text( colour="black", size = 25),
          axis.title.y = element_text( colour="black", size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(face="bold", size = 10))+
    scale_fill_manual(values = c("#BE433D","#459DA5"))+
    scale_color_manual(values =  c("#BE433D","#459DA5"))
  
  p.os + scale_fill_manual(values = c("white", "white"))
  
  p.save = p.os + scale_fill_manual(values = c("white", "white")) + 
    geom_step(aes(color = group), size =5, alpha = 0.8) + 
    geom_point(size = 5, shape = "+") +
    scale_color_manual(values =  c("#BE433D","#459DA5")) +
    annotate(geom = "text", label = paste0("p =", round(p.val, 15)), x = 1000, y = 0.7, size = 5) +
    annotate(geom = "text", label = paste0("p.sig =", round(p.val.signature, 15)), x = 1000, y = 0.75, size = 5) +
    ggtitle(label = seurat.obj.name, subtitle = paste("n=", nrow(survival.data))) +
    scale_y_continuous(limits = c(lower.limit, 1))
  
  
  return(p.save)
  
  
  
}




#final gene set 2
plot.list = list()
for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(score.name = "BMP17", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k])
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("BMP17.AUC.thirds", ".pdf"), width = 20, height = 12, 
         path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/5_BMP17_VST_AUC")
  
}

#BMPsurface
for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(score.name = "BMP.surface.AUC", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k])
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("BMP.surface.AUC", ".pdf"), width = 20, height = 12, 
         path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/5_BMP17_VST_AUC")
  
}


#TCF7LEF1
for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(score.name = "TCF7LEF1.AUC", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k], lower.limit = 0.4)
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("TCF7LEF1.AUC", ".pdf"), width = 20, height = 12, 
         path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/5_BMP17_VST_AUC")
  
}





#BMP surface

plot.list = list()
for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(score.name = "BMP17", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k])
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("BMP17", ".pdf"), width = 20, height = 12, 
         path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/5_BMP17_VST_AUC")
  
}






#import gene sig

consensus.bmp.sig.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/consensus.bmp.sig.nofilter.rds")
consensus.t.spec.sig.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/consensus.t.spec.sig.nofilter.rds")




#bmp consesnsus
ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9 <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9.rds")
ETP.Non.ETP.consensus.T.spec.sig.mean.log2FC.0.9 <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/ETP.Non.ETP.consensus.T.spec.sig.mean.log2FC.0.9.rds")

#write_csv(ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9, file = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9.csv")


#bmp consensus with scATAC merged TFs
bmp.scatac.tfs = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure4_Progenitor_Cell_Population/4F_intersect_DE/top.tfs.de.and.da.rds")
bmp.tfs.add = bmp.scatac.tfs %>% filter(cluster1 == "BMP-like-NR")
bmp.tfs.add %>% mutate(total.log10.q = log(p_val_adj, base = 10) + log(pv_adjust, base = 10)) %>% relocate(total.log10.q)

t.spec.tfs.add = bmp.scatac.tfs %>% filter(cluster1 != "BMP-like-NR")


#make list of sigs and seurats

p.features = list(consensus.bmp.sig.nofilter$gene, 
                  ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9$gene, 
                  c(consensus.bmp.sig.nofilter$gene, bmp.tfs.add$feature), 
                  c(ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9$gene, bmp.tfs.add$feature))

#names(p.features) = c("consensus.bmp.sig.nofilter (n=26)", "consensus.bmp.filter (n=17)", "consensus.bmp.sig.no.filter.plus.scATAC.tfs (n=38)", "consensus.bmp.filter.plus.scATAC.tfs (n=29)")
names(p.features) = c("consensus.bmp.sig.nofilter.26", "consensus.bmp.filter.17", "consensus.bmp.sig.no.filter.plus.scATAC.tfs.38", "consensus.bmp.filter.plus.scATAC.tfs.29")

n.features = list(consensus.t.spec.sig.nofilter$gene, 
                  ETP.Non.ETP.consensus.T.spec.sig.mean.log2FC.0.9$gene, 
                  c(consensus.t.spec.sig.nofilter$gene, t.spec.tfs.add$feature), 
                  c(ETP.Non.ETP.consensus.T.spec.sig.mean.log2FC.0.9$gene, t.spec.tfs.add$feature))
lapply(n.features, length)
#names(n.features) = c("consensus.t.spec.sig.nofilter (n=21)", "consensus.t.spec.filter (n=10)", "consensus.t.spec.sig.no.filter.plus.scATAC.tfs (n=28)", "consensus.t.spec.filter.plus.scATAC.tfs (n=17)")
names(n.features) = c("consensus.t.spec.sig.nofilter.21", "consensus.t.spec.filter.10", "consensus.t.spec.sig.no.filter.plus.scATAC.tfs.28", "consensus.t.spec.filter.plus.scATAC.tfs.17")



#BMP119





bulk.rna.object.umap20@meta.data %>% colnames()
bulk.rna.object.umap20@meta.data  = bulk.rna.object.umap20@meta.data %>% mutate(EFS = Event.free.Survival.in.Days,
                                                                                EFS.Censored = Censored.Event.free.Survival.in.Days,
                                                                                OS = Overall.Survival.in.Days,
                                                                                OS.Censored = Censored.Overall.Survival.in.Days,
                                                                                DFS = Disease.free.Survival.in.Days,
                                                                                DFS.Censored = Censored.Disease.free.Survival.in.Days,
                                                                                D29.MRD = Day.29.MRD %>% as.numeric(),
                                                                                D29.MRD.pos = Day.29.MRD > 0.01,
                                                                                EOC.MRD = End.of.Consolidation.MRD,
                                                                                D29.BM = Day.29.morphologic.Response, 
                                                                                cns.status = CNS.Status,
                                                                                event = Event.Type,
                                                                                ETP = ETP.STATUS, 
                                                                                dx.age = Age.at.Diagnosis.in.Days,
                                                                                dx.wbc = WBC.at.Diagnosis)

bulk.rna.object.umap20@meta.data$D29.MRD %>% class()
#0434 All
bulk.rna.object.umap20@meta.data = bulk.rna.object.umap20@meta.data %>% survival.conversion.for.plotting()

#saveRDS(bulk.rna.object.umap20, "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/1_harmonizemtx_make_seurat/vst.no.relapse.bulk.rna.umap.20.rds")

#0434 ETP
etp.bulk = subset(bulk.rna.object.umap20, ETP == "ETP")

#0434 Near-ETP
near.etp.bulk = subset(bulk.rna.object.umap20, ETP == "Near-ETP")

#0434 Non-ETP
non.etp.bulk = subset(bulk.rna.object.umap20, ETP == "Non-ETP")


bulk.rna.objects = list(bulk.rna.object.umap20, etp.bulk, near.etp.bulk, non.etp.bulk)
names(bulk.rna.objects) = c("0434 All", "0434 ETP", "0434 Near-ETP", "0434 Non-ETP")

#for each signature, stratify all 6 seurat objects

bulk.rna.object.umap20$ETP %>% table()



do.survival.given.gene.set = function(pos.features, neg.features = NA, score.name, seurat.object, seurat.obj.name, lower.limit = 0.6){
  
  print('scoring object')
  if(is.na(neg.features)){
    seurat.object <- signatureScore_zscore(seurat.object, 
                                           pfeatures = pos.features,
                                           score.name = score.name)
  }else{
    seurat.object <- signatureScore_zscore(seurat.object, 
                                           pfeatures = pos.features,
                                           nfeatures = neg.features,
                                           score.name = score.name)
  }
  
  
  print('fixing surv metadata')
  
  #seurat.object@meta.data = seurat.object@meta.data %>% survival.conversion.for.plotting()
  seurat.object@meta.data = seurat.object@meta.data %>% mutate(status.OS.plot = case_when(status.OS == 0 ~ "Alive",
                                                                                          status.OS == 1 ~ "Dead"))
  
  
  
  #survival plot
  library(survival)
  library(ggfortify)
  
  
  
  print('running cox')
  survival.data = seurat.object@meta.data
  
  
  survival.data$mpp.rank =  survival.data[[score.name]] %>% rank()
  survival.data$status.OS.plot %>% table()
  
  
  survival.data = survival.data %>% mutate(high.mpp = case_when(mpp.rank > (nrow(survival.data)*1/2) ~ "BMP High",
                                                                mpp.rank <= (nrow(survival.data)*1/2) ~ "BMP Low"))
  
  
  table(survival.data$high.mpp)
  survival.data$status.OS = as.numeric(survival.data$status.OS)
  os_km_fit <- survfit(Surv(time.OS, status.OS) ~ high.mpp , data=survival.data)
  cox.test = coxph(Surv(time.OS, status.OS) ~ high.mpp + D29.MRD, data=survival.data)
  cox.coefs.gene.and.mrd = summary(cox.test)
  
  if(nrow(survival.data) >100 ){
    p.val = cox.coefs.gene.and.mrd$logtest[3] %>% print()
  }else{
    p.val = cox.coefs.gene.and.mrd$sctest[3] %>% print()
    
  }
  
  
  
  if("cns.status" %in% colnames(survival.data)){
    cox.test2 = coxph(Surv(time.OS, status.OS) ~ survival.data[[score.name]] + D29.MRD + cns.status, data=survival.data)
    cox.test.summary = summary(cox.test2)
    p.val.signature = cox.test.summary$coefficients[1,5] %>% print()
    
  }else{
    cox.test2 = coxph(Surv(time.OS, status.OS) ~ survival.data[[score.name]] + D29.MRD, data=survival.data)
    cox.test.summary = summary(cox.test2)
    p.val.signature = cox.test.summary$coefficients[1,5] %>% print()
    
  }
  
  
  #ccp.val = cox.coefs.gene.and.mrd$sctest[3] %>% print()
  
  p.os = autoplot(os_km_fit)+ 
    labs(x = "OS (Days) ", y = "Survival Probability", 
         title = "") + 
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_text( colour="black", size = 25),
          axis.title.y = element_text( colour="black", size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(face="bold", size = 10))+
    scale_fill_manual(values = c("#BE433D","#459DA5"))+
    scale_color_manual(values =  c("#BE433D","#459DA5"))
  
  p.os + scale_fill_manual(values = c("white", "white"))
  
  p.save = p.os + scale_fill_manual(values = c("white", "white")) + 
    geom_step(aes(color = group), size =5, alpha = 0.8) + 
    geom_point(size = 5, shape = "+") +
    scale_color_manual(values =  c("#BE433D","#459DA5")) +
    annotate(geom = "text", label = paste0("p =", round(p.val, 15)), x = 1000, y = 0.7, size = 5) +
    annotate(geom = "text", label = paste0("p.sig =", round(p.val.signature, 15)), x = 1000, y = 0.75, size = 5) +
    ggtitle(label = seurat.obj.name, subtitle = paste("n=", nrow(survival.data))) +
    scale_y_continuous(limits = c(lower.limit, 1))
  
  
  return(p.save)
  
  
  
}




#final gene set 2
plot.list = list()
for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(pos.features =  p.features$consensus.bmp.filter.17,
                                              score.name = "consensus.bmp.17", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k])
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("consensus.bmp.17", ".pdf"), width = 20, height = 12, path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/2_do_BMP_17_survival")
  
}

#final gene set 2
plot.list = list()
for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(pos.features =  p.features$consensus.bmp.filter.17,
                                              neg.features = n.features$consensus.t.spec.filter.10,
                                              score.name = "consensus.bmp.17.t.spec.sig.filter.10", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k])
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("consensus.bmp.17.t.spec.sig.filter.10", ".pdf"), width = 20, height = 12, path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/2_do_BMP_17_survival")
  
}

consensus.bmp.surface.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/BMP_surface_phenotype/consensus.bmp.surface.nofilter.rds")
consensus.t.spec.surface.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/BMP_surface_phenotype/consensus.t.spec.surface.nofilter.rds")

#final gene set 2
plot.list = list()
for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(pos.features =  consensus.bmp.surface.nofilter$coding,
                                              neg.features = consensus.t.spec.surface.nofilter$coding,
                                              score.name = "BMP.surface.9", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k])
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("BMP.surface.9", ".pdf"), width = 15, height = 12, path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/2_do_BMP_17_survival")
  
}

BMP.119.ETP <- read_csv("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/Figure2/119_gene_deg.sig.csv")

for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(pos.features =  BMP.119.ETP$gene[BMP.119.ETP$enriched == "BMP-like-NR"],
                                              neg.features = BMP.119.ETP$gene[BMP.119.ETP$enriched == "T-specified-R"],
                                              score.name = "BMP.119.ETP", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k], lower.limit = 0)
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("BMP.119.ETP", ".pdf"), width = 15, height = 12, path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/2_do_BMP_17_survival")
  
}





BMP.119.ETP <- read_csv("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/Figure2/119_gene_deg.sig.csv")

for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(pos.features =  BMP.119.ETP$gene[BMP.119.ETP$enriched == "BMP-like-NR"],
                                              neg.features = BMP.119.ETP$gene[BMP.119.ETP$enriched == "T-specified-R"],
                                              score.name = "BMP.119.ETP", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k], lower.limit = 0)
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("BMP.119.ETP", ".pdf"), width = 15, height = 12, path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/2_do_BMP_17_survival")
  
}

#TCF7-LEF1



t.all.40.eps=read_tsv("/mnt/isilon/tan_lab/yuw1/R_work_dir/ETP_ALL/EP_Prediction2/regrRes4_EP_overall.txt")
edges = read_tsv("/mnt/isilon/tan_lab/yuw1/R_work_dir/ETP_ALL/EP_Prediction2/TRN_edges_enrichedTFs_in_Near-ETP_stringentPlus.txt")
nodes = read_tsv("/mnt/isilon/tan_lab/yuw1/R_work_dir/ETP_ALL/EP_Prediction2/TRN_nodes_enrichedTFs_in_Near-ETP_stringentPlus.txt")
ep.predictions.t.all.40 <- readRDS("/mnt/isilon/tan_lab/yuw1/R_work_dir/ETP_ALL/EP_Prediction2/regrRes4ep_prediction.rds")

t.all.40.eps %>% filter(gene_name == "CD5") %>% arrange(fdr)

edges %>% filter(gene_name == "CD5")


tcf7.lef1.targets = edges %>% filter(TF %in% c("TCF7", "LEF1"))
tcf7.lef1.targets.up = tcf7.lef1.targets %>% filter(direction == "up")
tcf7.lef1.targets.down = tcf7.lef1.targets %>% filter(direction == "down")

tcf7.lef1.cd5.activtors = edges %>% filter(gene_name %in% c("TCF7", "LEF1", "CD5")) 
tcf7.lef1.activtors = edges %>% filter(gene_name %in% c("TCF7", "LEF1")) 


for(k in 1:length(bulk.rna.objects)){
  print(k)
  
  plot.list[[k]] = do.survival.given.gene.set(pos.features =  c(tcf7.lef1.targets.up$gene_name, tcf7.lef1.cd5.activtors$TF) %>% unique(),
                                              #neg.features = tcf7.lef1.targets.down$gene_name %>% unique(),
                                              score.name = "tcf7lef1up", 
                                              seurat.object = bulk.rna.objects[[k]], 
                                              seurat.obj.name = names(bulk.rna.objects)[k])
  
  
  to.save = cowplot::plot_grid(plotlist = plot.list)
  ggsave(to.save, filename = paste0("tcf7lef1up-down-with-activators", ".pdf"), width = 15, height = 12, 
         path = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/2_do_BMP_17_survival")
  
}

