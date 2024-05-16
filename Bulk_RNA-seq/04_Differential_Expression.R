library(limma)

setwd("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/9_remake_NonETP_MRD+_vsMRD-_volcano")

load("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/data_sharing/TALL_X01_normalized_1362.Rdata")

eset_counts.old <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/bulkRNA/13_FinalMtx/eset_counts.rds")

rownames(patient_anno) = patient_anno$USI

patient_anno  = patient_anno %>% mutate(EFS = Event.free.Survival.in.Days,
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

rownames(gene_anno) = gene_anno$geneID
rownames(gene_anno)

eset.new <- ExpressionSet(assayData =  unfiltered.counts, #eset object
                      phenoData = AnnotatedDataFrame(patient_anno),
                      featureData = AnnotatedDataFrame(gene_anno))






eset = eset.new

eset <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/9_remake_NonETP_MRD+_vsMRD-_volcano/1_02_eset_counts.rds")

#do deg testing

pData(eset) = pData(eset) %>% mutate(response= case_when(D29.MRD >= 1 ~ "mrd.over.1",
                                                         D29.MRD >= 0.01 & D29.MRD < 1 ~ "mrd.01.1",
                                                         D29.MRD < 0.01 ~ "no.mrd"))

#first try binary test
pData(eset) = pData(eset) %>% mutate(response= case_when(D29.MRD >= 0.01 ~ "mrd.positive",
                                                         D29.MRD < 0.01 ~ "no.mrd"))
etp.eset = eset[,eset$ETP == "ETP"]

#start DE testing----
group <- with(pData(eset), paste(ETP, response, sep = ".")) #cat-ing groups for design matrix
group = gsub(pattern = "-", replacement = ".", x = group) #replace -

group <- factor(group) %>% print()
table(group)


design <- model.matrix(~0 + group) %>% print() #design matrix
colnames(design) <- levels(group)
colSums(design) #checking design matrix



#higher.filter.threshold = rowSums(exprs(eset) > 25) > 50
#higher.filter.threshold = rowSums(log2.cpm.filtered.norm.df > 2) > 20

#v.eset <- voom(etp.eset, design, plot = TRUE) #running voom on eset --> log2cpm counts and weights.


#keepers
unfiltered.data.mtx = eset@assayData$exprs %>% print()
table(rowSums(unfiltered.data.mtx>1)>=50)

keepers <- rowSums(unfiltered.data.mtx>1)>=50
table(keepers)


v.eset <- voom(eset[keepers,], design, plot = TRUE) #running voom on eset --> log2cpm counts and weights.



cm2 <- makeContrasts(ETP.MRD.Over.0.01.vs.No.MRD = ETP.mrd.positive - ETP.no.mrd,
                     NearETP.MRD.Over.0.01.vs.No.MRD = Near.ETP.mrd.positive - Near.ETP.no.mrd,
                     NonETP.MRD.Over.0.01.vs.No.MRD = Non.ETP.mrd.positive - Non.ETP.no.mrd,
                     levels = design)

fit <- lmFit(v.eset, design)
fit_all <- contrasts.fit(fit, contrasts = cm2)
fit_all <- eBayes(fit_all)

summary(decideTests(fit_all, adjust.method="BH", p.value=0.05))

results_all_p025 <- decideTests(fit_all, adjust.method="BH", p.value=0.25, method = "separate")

#2 or more, same direction, down
dges <- results_all_p025 %>% as.data.frame.matrix() %>% 
  as.data.frame(rownames = "geneID")


topTable(fit_all, adjust ="BH", coef=1) %>% arrange(adj.P.Val) %>% head(50)

setwd("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/9_remake_NonETP_MRD+_vsMRD-_volcano")
figure.directory = "/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/9_remake_NonETP_MRD+_vsMRD-_volcano"
#make volcanoplots----
volcano.fig.dir = paste0(figure.directory, "/VolcanoPlots_mrdPos_vs_mrdNeg") %>% print()
dir.create(volcano.fig.dir)

limma.list = list()
library(EnhancedVolcano)

for(i in 1:length(colnames(fit_all))){
  
  comparison = colnames(fit_all)[i] %>% print()
  
  
  limma.list[[i]] = topTable(fit_all, coef = i, adjust ="BH", p.value = 0.5, number = 100000) %>%
    as_tibble() %>% print()
  
  if(nrow(limma.list[[i]]) ==0){
    print("no sig genes")
    
  }else{
    
    #volcanoplot(fit_all, coef = i, style = "p-value", highlight = 20, names = fit$genes$geneID, hl.col="blue",
    #xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)
    #ggsave(path = volcano.fig.dir, paste0("limma_default_", comparison, ".png"), height = 10, width = 10)
    
    
    
    
    top.genes = limma.list[[i]]$geneID[limma.list[[i]]$logFC>0] %>% head(25)
    top.genes = c(top.genes, limma.list[[i]]$geneID[limma.list[[i]]$logFC<0] %>% head(25))
    
    ylim = limma.list[[i]]$adj.P.Val %>% min()
    
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneSymbol,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = comparison,
                    #selectLab = top.genes,
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    pointSize = 3.0,
                    #ylim = c(0, -log10(10e-5)),
                    ylim = c(0, -log10(ylim / 1.1)),
                    labSize = 3, drawConnectors = T) 
    
    ggsave(path = volcano.fig.dir, paste0(comparison, "_ALL.pdf"), height = 10, width = 10)
    
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneSymbol,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = comparison,
                    #selectLab = top.genes,
                    pCutoff = 0.05,
                    FCcutoff = 0.5,
                    pointSize = 3.0,
                    labSize = 3, drawConnectors = F, 
                    ylim = c(0, -log10(ylim / 1.1))) 
    
    ggsave(path = volcano.fig.dir, paste0(comparison, "_Default.pdf"), height = 10, width = 10)
    
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneSymbol,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = paste(comparison, "0.5FC & 0.25FDR"),
                    #selectLab = top.genes,
                    pCutoff = 0.25,
                    FCcutoff = 0.5,
                    pointSize = 3.0,
                    labSize = 3, drawConnectors = F, 
                    ylim = c(0, -log10(ylim / 1.1))) 
    
    ggsave(path = volcano.fig.dir, paste0(comparison, "_0.5FC-0.25Pval.pdf"), height = 10, width = 10)
    
    
  }
  names(limma.list)[[i]] = comparison
  
}


#nonetp
source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/JX_scFunctions.R')
binary.comparisons = limma.list #save limma list
human.tfs.and.cisbp = get.human.tfs()
gene.de = binary.comparisons$NonETP.MRD.Over.0.01.vs.No.MRD
gene.de$comparison = "NonETP.MRD.Over.0.01.vs.No.MRD"
gene.de$enriched = ifelse(gene.de$logFC > 0, yes = "MRD+", no = "MRD-")
gene.de$tf = gene.de$geneSymbol %in% human.tfs.and.cisbp$human.tfs.and.cisbp

genes.only = gene.de %>% filter(tf == F) %>% left_join(eset@featureData@data)
bulk.volcanos = list()
bulk.volcanos$gene = ggplot(genes.only %>% filter(abs(logFC) > 1, adj.P.Val < 1e-10, bioType %in% c("protein_coding", "TR_C_gene")), 
                            aes(x = -logFC, y = -log(adj.P.Val, base = 10),
                                color = enriched,
                                label = geneSymbol, 
                                shape = tf)) + 
  geom_point(data = genes.only %>% filter(bioType == "protein_coding"), color = "light gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 50) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  #scale_color_manual(values = c("#BE433D","#DDA164","#B89B74")) +
  scale_color_manual(values =  c("black","dark red"))+
  #geom_vline(xintercept = 0.5, linetype = "dashed")+
  #geom_vline(xintercept = -0.5, linetype = "dashed")+
  #geom_hline(yintercept = -log(1e-25, base = 10), linetype = "dashed") + 
  facet_wrap(~comparison, nrow = 1) + 
  theme(legend.position = "none")

tf.de = gene.de %>% filter(tf == T)
bulk.volcanos$tf = ggplot(tf.de %>% filter(abs(logFC) > 0.25, adj.P.Val < 1e-05), 
                          aes(x = -logFC, y = -log(adj.P.Val, base = 10),
                              color = enriched,
                              label = geneSymbol, 
                              shape = tf)) + 
  geom_point(data = tf.de, color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 50) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  #scale_color_manual(values = c("#BE433D","#DDA164","#B89B74")) +
  scale_color_manual(values =  c("black","dark red"))+
  #geom_vline(xintercept = 0.5, linetype = "dashed")+
  #geom_vline(xintercept = -0.5, linetype = "dashed")+
  #geom_hline(yintercept = -log(1e-25, base = 10), linetype = "dashed") + 
  facet_wrap(~comparison, nrow = 1) + 
  theme(legend.position = "none")





ggplot(gene.de %>% filter(abs(logFC) > 1, adj.P.Val < 1e-10, bioType %in% c("protein_coding", "TR_C_gene")), 
       aes(x = -logFC, y = -log(adj.P.Val, base = 10),
           color = enriched,
           label = geneSymbol, 
           shape = tf)) + 
  geom_point(data = gene.de %>% filter(bioType == "protein_coding"), color = "light gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 50) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  #scale_color_manual(values = c("#BE433D","#DDA164","#B89B74")) +
  scale_color_manual(values =  c("black","dark red"))+
  #geom_vline(xintercept = 0.5, linetype = "dashed")+
  #geom_vline(xintercept = -0.5, linetype = "dashed")+
  #geom_hline(yintercept = -log(1e-25, base = 10), linetype = "dashed") + 
  facet_wrap(~comparison, nrow = 1) + 
  theme(legend.position = "none")


#
#make DEG by three gruops----
pData(eset) = pData(eset) %>% mutate(response= case_when(D29.MRD >= 1 ~ "mrd.over.1",
                                                         D29.MRD >= 0.01 & D29.MRD < 1 ~ "mrd.01.1",
                                                         D29.MRD < 0.01 ~ "no.mrd"))

#etp.eset = eset[,eset$ETP == "ETP"]

#start DE testing
group <- with(pData(eset), paste(ETP, response, sep = ".")) #cat-ing groups for design matrix
group = gsub(pattern = "-", replacement = ".", x = group) #replace -

group <- factor(group) %>% print()
table(group)


design <- model.matrix(~0 + group) %>% print() #design matrix
colnames(design) <- levels(group)
colSums(design) #checking design matrix



#higher.filter.threshold = rowSums(exprs(eset) > 25) > 50
#higher.filter.threshold = rowSums(log2.cpm.filtered.norm.df > 2) > 20

#v.eset <- voom(etp.eset, design, plot = TRUE) #running voom on eset --> log2cpm counts and weights.

#v.eset <- voom(eset[higher.filter.threshold,], design, plot = TRUE) #running voom on eset --> log2cpm counts and weights.


cm <- makeContrasts(ETP.MRD.Over.1.vs.No.MRD =  ETP.mrd.over.1 - ETP.no.mrd, #first Combo vs DMSO
                    ETP.MRD.0.01.to.1.vs.No.MRD = ETP.mrd.01.1 - ETP.no.mrd,
                    ETP.MRD.Over.1.vs.MRD.0.01.to.1=  ETP.mrd.over.1- ETP.mrd.01.1 ,
                    
                    Near.ETP.MRD.Over.1.vs.No.MRD =  Near.ETP.mrd.over.1 - Near.ETP.no.mrd, #first Combo vs DMSO
                    Near.ETP.MRD.0.01.to.1.vs.No.MRD = Near.ETP.mrd.01.1 - Near.ETP.no.mrd,
                    Near.ETP.MRD.Over.1.vs.MRD.0.01.to.1=  Near.ETP.mrd.over.1- Near.ETP.mrd.01.1 ,
                    
                    Non.ETP.MRD.Over.1.vs.No.MRD =  Non.ETP.mrd.over.1 - Non.ETP.no.mrd, #first Combo vs DMSO
                    Non.ETP.MRD.0.01.to.1.vs.No.MRD = Non.ETP.mrd.01.1 - Non.ETP.no.mrd,
                    Non.ETP.MRD.Over.1.vs.MRD.0.01.to.1=  Non.ETP.mrd.over.1- Non.ETP.mrd.01.1 ,
                    levels = design)

fit <- lmFit(v.eset, design)
fit_all <- contrasts.fit(fit, contrasts = cm)
fit_all <- eBayes(fit_all)

summary(decideTests(fit_all, adjust.method="BH", p.value=0.05))




#make volcanoplots----
volcano.fig.dir = paste0(figure.directory, "/VolcanoPlots_Stratified/") #make a new one
dir.create(volcano.fig.dir)

limma.list = list()
#library(EnhancedVolcano)

for(i in 1:length(colnames(fit_all))){
  
  comparison = colnames(fit_all)[i] %>% print()
  
  
  limma.list[[i]] = topTable(fit_all, coef = i, number=300000, adjust ="BH", p.value = 0.5) %>%
    as_tibble() %>% print()
  
  if(nrow(limma.list[[i]]) ==0){
    print("no sig genes")
    
  }else{
    
    #volcanoplot(fit_all, coef = i, style = "p-value", highlight = 20, names = fit$genes$geneID, hl.col="blue",
    #xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)
    #ggsave(path = volcano.fig.dir, paste0("limma_default_", comparison, ".png"), height = 10, width = 10)
    
    
    
    
    top.genes = limma.list[[i]]$geneID[limma.list[[i]]$logFC>0] %>% head(25)
    top.genes = c(top.genes, limma.list[[i]]$geneID[limma.list[[i]]$logFC<0] %>% head(25))
    
    ylim = limma.list[[i]]$adj.P.Val %>% min()
    
    pdf(paste0(volcano.fig.dir, comparison, "_3Plots.pdf"), height = 10, width = 10)
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneSymbol,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = paste(comparison, "0.5FC & 0.25FDR"),
                    #selectLab = top.genes,
                    pCutoff = 0.25,
                    FCcutoff = 0.5,
                    pointSize = 3.0,
                    labSize = 3, drawConnectors = F, 
                    ylim = c(0, -log10(ylim / 1.1)))%>% print()
    
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneSymbol,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = paste(comparison, "ConnectingLabels"),
                    #selectLab = top.genes,
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    pointSize = 3.0,
                    #ylim = c(0, -log10(10e-5)),
                    ylim = c(0, -log10(ylim / 1.1)),
                    labSize = 3, drawConnectors = T) %>% print()
    
    
    
    EnhancedVolcano(limma.list[[i]],
                    lab = limma.list[[i]]$geneSymbol,
                    x = 'logFC',
                    y = 'adj.P.Val',
                    title = paste(comparison, "DefaultPlot"),
                    #selectLab = top.genes,
                    pCutoff = 0.05,
                    FCcutoff = 0.5,
                    pointSize = 3.0,
                    labSize = 3, drawConnectors = F, 
                    ylim = c(0, -log10(ylim / 1.1))) %>% print()
    
    
    dev.off()
    
    
  }
  names(limma.list)[[i]] = comparison
  
}


mrd.over.1.vs.mrd.neg = limma.list
human.tfs.and.cisbp = get.human.tfs()
gene.de2 = mrd.over.1.vs.mrd.neg$Non.ETP.MRD.Over.1.vs.No.MRD
gene.de2$comparison = "Non.ETP.MRD.Over.1.vs.No.MRD"
gene.de2$enriched = ifelse(gene.de2$logFC > 0, yes = "MRD>1", no = "MRD-")
gene.de2$tf = gene.de2$geneSymbol %in% human.tfs.and.cisbp$human.tfs.and.cisbp

genes.only2 = gene.de2 %>% filter(tf == F) %>% left_join(eset@featureData@data)
#bulk.volcanos = list()
bulk.volcanos$gene2 = ggplot(genes.only2 %>% filter(abs(logFC) > 1, adj.P.Val < 1e-10, bioType %in% c("protein_coding", "TR_C_gene")), 
                             aes(x = -logFC, y = -log(adj.P.Val, base = 10),
                                 color = enriched,
                                 label = geneSymbol, 
                                 shape = tf)) + 
  geom_point(data = genes.only2 %>% filter(bioType == "protein_coding"), color = "light gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 50) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  #scale_color_manual(values = c("#BE433D","#DDA164","#B89B74")) +
  scale_color_manual(values =  c("black","dark red"))+
  #geom_vline(xintercept = 0.5, linetype = "dashed")+
  #geom_vline(xintercept = -0.5, linetype = "dashed")+
  #geom_hline(yintercept = -log(1e-25, base = 10), linetype = "dashed") + 
  facet_wrap(~comparison, nrow = 1) + 
  theme(legend.position = "none")
bulk.volcanos$gene2
tf.de2 = gene.de2 %>% filter(tf == T)
bulk.volcanos$tf2 = ggplot(tf.de2 %>% filter(abs(logFC) > 0.25, adj.P.Val < 1e-05), 
                           aes(x = -logFC, y = -log(adj.P.Val, base = 10),
                               color = enriched,
                               label = geneSymbol, 
                               shape = tf)) + 
  geom_point(data = tf.de2, color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 50) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)")+
  #scale_color_manual(values = c("#BE433D","#DDA164","#B89B74")) +
  scale_color_manual(values =  c("black","dark red"))+
  #geom_vline(xintercept = 0.5, linetype = "dashed")+
  #geom_vline(xintercept = -0.5, linetype = "dashed")+
  #geom_hline(yintercept = -log(1e-25, base = 10), linetype = "dashed") + 
  facet_wrap(~comparison, nrow = 1) + 
  theme(legend.position = "none")
bulk.volcanos$tf2

getwd()
pdf("bulk.de.plots.pdf", height = 25, width = 25)
bulk.volcanos
dev.off()

