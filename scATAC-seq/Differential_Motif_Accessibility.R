source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/JX_scFunctions.R')
setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL")

## prepare gene expression fraction in the corresponding cell populations ####
seurat.rna = readRDS("SCENICPlus/etp.25.g1.downsample.1711.each.patient.rds")
table(seurat.rna$comparison.bmp.vs.t.specified)

blasts.seurat.ann.metadata = seurat.rna@meta.data %>% 
  mutate(comparison.bmp.vs.t.specified = case_when(response.binary == "mrd.over.01" & predicted.cell.type.short %in% c("HSPC", "LMPP", "CLP", "ETP") & ETP == "ETP" ~ "BMP-like-NR",
                                                   response.binary == "no.mrd" & predicted.cell.type.short %in% c("Pro-T", "Pre-T") & ETP == "ETP" ~ "T-specified-R", 
                                                   response.binary == "no.mrd" & predicted.cell.type.short %in% c("HSPC", "LMPP", "CLP", "ETP") & ETP == "ETP" ~ "BMP-like-R",
                                                   response.binary == "mrd.over.01" & predicted.cell.type.short %in% c("Pro-T", "Pre-T") & ETP == "ETP" ~ "T-specified-NR"))
seurat.rna = AddMetaData(seurat.rna, metadata = blasts.seurat.ann.metadata)
table(seurat.rna$comparison.bmp.vs.t.specified)

Ctypes = seurat.rna$comparison.bmp.vs.t.specified
Ctypes[is.na(Ctypes)] = "Other"
mtx = seurat.rna@assays$RNA@counts

#load chromVAR results
jx.save = '/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/18_chromVAR_etp-near-non/'
load(paste0(jx.save, 'chromVARObj/zscore_dscore_allPeaks_ETP.RData')) #loads chromVAR results

#get annotations of each cell
seurat.atac <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds')
seurat.atac@meta.data = seurat.atac@meta.data %>% mutate(sample.name.cell.barcode = paste0(sample.name, "-", cell.barcode))

source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/JX_scFunctions.R')
setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL")

all.etp.60.mdata = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/02_projectv4_atac/2022-03-24_ProjectionFiles_10_PCs_min_dist0.5_7171_VAPs/umap_coords_2_with_metadata_and_maturity.RDS")

ann.atac = data.table(seurat.atac@meta.data, keep.rownames = T)

#add new metadata
all.etp.60.mdata = all.etp.60.mdata %>% filter(cell.name %in% seurat.atac$sample.name.cell.barcode) %>% left_join(seurat.atac@meta.data %>% select("chromVAR.cell.name" = cell.name, 
                                                                                                                                                   "cell.name" = sample.name.cell.barcode,
                                                                                                                                                   response.binary, ETP))
ann.atac = all.etp.60.mdata
ann.atac = ann.atac %>% 
  mutate(comparison.bmp.vs.t.specified = case_when(response.binary == "mrd.over.01" & predicted.cell.type.short %in% c("HSPC", "LMPP", "CLP", "ETP") & ETP == "ETP" ~ "BMP-like-NR",
                                                   response.binary == "no.mrd" & predicted.cell.type.short %in% c("Pro-T", "Pre-T") & ETP == "ETP" ~ "T-specified-R", 
                                                   response.binary == "no.mrd" & predicted.cell.type.short %in% c("HSPC", "LMPP", "CLP", "ETP") & ETP == "ETP" ~ "BMP-like-R",
                                                   response.binary == "mrd.over.01" & predicted.cell.type.short %in% c("Pro-T", "Pre-T") & ETP == "ETP" ~ "T-specified-NR"))
ann.atac$comparison.bmp.vs.t.specified[is.na(ann.atac$comparison.bmp.vs.t.specified)] = "Other"

ann.atac$comparison.bmp.vs.t.specified %>% table()

ann.atac %>% group_by(ETP, response.binary) %>% summarize(n = length(ETP), 
                                                          bmp = sum(predicted.cell.type.short %in% c("HSPC", "LMPP", "CLP", "ETP")), 
                                                          t.spec = sum(predicted.cell.type.short %in% c("Pro-T", "Pre-T")),
                                                          t.lineage = sum(predicted.cell.type.short %in% c("alpha-beta", "alpha-beta(mature)", "Naive_T", "DP", "EffectorT")))


cls.compare = data.table(barcode = ann.atac$chromVAR.cell.name,
                         cluster = ann.atac$comparison.bmp.vs.t.specified) %>% arrange(cluster)
cls.compare$cluster %>% table()


###
#1
comparisons = c("T-specified-R", "T-specified-NR")

#2
comparisons = c("BMP-like-R", "BMP-like-NR") 

#3
comparisons = c("T-specified-R", "BMP-like-NR")

#4
comparisons = c("BMP-like-R", "T-specified-R")

#5
comparisons = c("BMP-like-R", "T-specified-NR")
####

efreq.1 <- rowMeans(mtx[, Ctypes == comparisons[2]] > 0)
efreq.2 <- rowMeans(mtx[, Ctypes == comparisons[1]] > 0)


#make tibbles of expression
efreq.1.tibble = tibble(feature = names(efreq.1),
                        expression = efreq.1)

efreq.2.tibble = tibble(feature = names(efreq.2),
                        expression = efreq.2)


#load function
source('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/17_run_chromVAR_r_vs_nr/scDataAnalysis_Utilities_simp.R')
setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL")

#load functions
runDiffMotifEnrich.with.median = function(mtx_score, clusters, test = 'wilcox',
                                          fdr = 0.01, topn = 20,
                                          min_frac_per_cluster = 0.2,
                                          max_cell_per_clust = 300,
                                          test.mode = 'one-rest', control = NULL){
  set.seed(2020)
  clusters$cluster = as.character(clusters$cluster)
  cls = unique(clusters$cluster)
  res = NULL
  features = rownames(mtx_score)
  mtx_score = mtx_score[, clusters$barcode]
  for(cluster0 in cls){
    bc0 = clusters[cluster == cluster0]$barcode
    mtx1 = mtx_score[, colnames(mtx_score) %in% bc0]
    if(test.mode == 'one-rest') {
      mtx2 = mtx_score[, !colnames(mtx_score) %in% bc0]
    }
    if(test.mode == 'control'){
      if(is.null(control)) stop('Should specific control cluster!')
      bc2 = clusters[cluster %in% control]$barcode
      mtx2 = mtx_score[, colnames(mtx_score) %in% bc2]
      if(cluster0 %in% control) next
    }
    mu1 = sapply(1:length(features), function(x) mean(mtx1[x, ]))
    mu2 = sapply(1:length(features), function(x) mean(mtx2[x, ]))
    
    median1 = sapply(1:length(features), function(x) median(mtx1[x, ]))
    median2 = sapply(1:length(features), function(x) median(mtx2[x, ]))
    
    
    pvs = rep(0.5, length(features))
    
    for(x in 1:length(features)){
      a1 = mtx1[x, ]
      a2 = mtx2[x, ]
      if(mu1[x] <= 0) next
      #if(mu1[x] <= mu2[x]) next
      if(length(which(!is.na(a1))) < 2 || length(which(!is.na(a2))) < 2) next
      if(mean(a1 > 0) < min_frac_per_cluster) next
      if(!is.null(max_cell_per_clust)){
        if(length(a1) > max_cell_per_clust) a1 = sample(a1, max_cell_per_clust)
        if(length(a2) > max_cell_per_clust) a2 = sample(a2, max_cell_per_clust)
      }
      
      if(test == 'wilcox') pvs[x] = wilcox.test(a1, a2, alternative = 'greater')$p.value
      if(test == 't') pvs[x] = t.test(a1, a2, alternative = 'greater')$p.value
      
    }
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster1' = cluster0,
                      'mean1' = mu1, 'mean0' = mu2, 
                      'median1' = median1, median2 = median2,
                      'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    res0 = res0[mean1 > 0]
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
} #for running DA



align.mean.median.chromVAR.da = function(res.diff.raw, exp.threshold = 0.2){
  #res.diff.raw = rbind(res.diff.raw.top.1000.with.median.etp.only, res.diff.raw.top.1000.with.median.near.non)
  res.diff.raw[, 'delta.mean' := mean1 - mean0]
  res.diff.raw[, 'delta.median' := median1 - median2]
  
  res.diff.raw = res.diff.raw[order(-delta.median)]
  #res.diff$cluster1 %>% table()
  
  #these lines are for filtering by expression. hash out if want to do later.
  #res.diff.raw.etp = res.diff.raw[cluster1 == 'ETP' & feature %in% names(which(efreq.etp >= exp.threshold))]
  #res.diff.raw.near.etp = res.diff.raw[cluster1 == 'Near-ETP' & feature %in% names(which(efreq.near.etp >= exp.threshold))]
  #res.diff.raw.non.etp = res.diff.raw[cluster1 == 'Non-ETP' & feature %in% names(which(efreq.non.etp >= exp.threshold))]
  #res.diff.raw.near.non.etp = res.diff.raw[cluster1 == 'Near-Non-ETP' & feature %in% names(which(efreq.near.non.etp >= exp.threshold))]
  
  
  res.diff.raw.two = res.diff.raw[cluster1 == comparisons[2]]
  res.diff.raw.one = res.diff.raw[cluster1 == comparisons[1]]
  
  
  print('adding expression')
  res.diff.raw.two=res.diff.raw.two %>% left_join(efreq.1.tibble)
  res.diff.raw.one= res.diff.raw.one %>% left_join(efreq.2.tibble)
  
  
  
  res.diff.raw.exp.filtered = rbind(res.diff.raw.two, res.diff.raw.one)
  res.diff.raw.exp.filtered
  
  res.diff.raw.exp.filtered.mean.median.alignment.filtered = res.diff.raw.exp.filtered %>% 
    mutate(mean.minus.median = delta.mean - delta.median,
           mean.div.median = delta.mean/delta.median) %>% 
    #filter(mean.div.median <= 1.15, mean.div.median >= 0.85 ) %>% 
    print(n=30)
  
  return(res.diff.raw.exp.filtered.mean.median.alignment.filtered)
} #for filtering by expression and alignment of mean and median

#perform DA

da.vs.one = runDiffMotifEnrich.with.median(dscore.all, 
                                           clusters = cls.compare, fdr = 0.01,
                                           max_cell_per_clust = 1500,
                                           topn = 1000, min_frac_per_cluster = 0.1,
                                           control = comparisons[1],
                                           test.mode = 'control')


da.vs.two = runDiffMotifEnrich.with.median(dscore.all, 
                                           clusters = cls.compare, fdr = 0.01,
                                           max_cell_per_clust = 1500,
                                           topn = 1000, min_frac_per_cluster = 0.1,
                                           control = comparisons[2],
                                           test.mode = 'control')

da.vs.one = da.vs.one %>% filter(cluster1 %in% c(comparisons[1], comparisons[2]))
da.vs.two = da.vs.two%>% filter(cluster1 %in% c(comparisons[1], comparisons[2]))

da.vs.one = da.vs.one %>% align.mean.median.chromVAR.da()
da.vs.two = da.vs.two %>% align.mean.median.chromVAR.da()


#make plot

da.two.vs.one = rbind(da.vs.one, da.vs.two)
da.two.vs.one$mean.div.median = da.two.vs.one$mean.div.median %>% round(2)
da.two.vs.one$cluster1 %>% table()

da.plots = list()
da.plots$da.two.vs.one.stringent = ggplot(da.two.vs.one %>% 
                                            filter(delta.median > 0.005, 
                                                   mean.div.median <= 1.1, 
                                                   mean.div.median >= 0.9, 
                                                   expression > 0.3,
                                                   pv_adjust < 1e-10), 
                                          aes(x = ifelse(cluster1 == comparisons[2], 
                                                         yes = -delta.median, no = delta.median), 
                                              y = -log(pv_adjust, base = 10),
                                              color = cluster1,
                                              label = feature)) + 
  geom_point(data = da.two.vs.one %>% filter(expression > 0.3,
                                             mean.div.median <= 1.1, 
                                             mean.div.median >= 0.9), color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 20) + 
  theme_bw() + 
  labs(x = "Accessibility Difference (delta chromVAR deviation)", 
       y = "-log(pv_adjust)")+
  scale_color_manual(values =  c("#BD20BF","#9DA52F")) +
  theme(legend.position = "none")



da.plots$da.two.vs.one.loose = ggplot(da.two.vs.one %>% 
                                        filter(delta.median > 0.0025, 
                                               mean.div.median < 1.3, 
                                               mean.div.median > 0.7, 
                                               expression > 0.2,
                                               pv_adjust < 1e-05), 
                                      aes(x = ifelse(cluster1 == comparisons[2], 
                                                     yes = -delta.median, no = delta.median), 
                                          y = -log(pv_adjust, base = 10),
                                          color = cluster1,
                                          label = feature)) + 
  geom_point(data = da.two.vs.one %>% filter(expression > 0.2,
                                             mean.div.median < 1.2, 
                                             mean.div.median > 0.8), color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 20) + 
  theme_bw() + 
  labs(x = "Accessibility Difference (delta chromVAR deviation)", 
       y = "-log(pv_adjust)")+
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  theme(legend.position = "none")

getwd()

#Saving plots 
#1
pdf("Responder_Figures/TSpecified_R_vs_NR_da.plots.downsample1500.pdf", height = 4, width = 4)
da.plots
dev.off()

#2
pdf("Responder_Figures/BMP-like_R_vs_BMP-like-NR_da.plots.downsample1500.pdf", height = 4, width = 4)
da.plots
dev.off()

#3
pdf("Responder_Figures/TSpecified_R_vs_BMP-like-NR_da.plots.downsample1500.pdf", height = 4, width = 4)
da.plots
dev.off()

#4
pdf("Responder_Figures/BMP-like_R_vs_TSpecified_R_da.plots.downsample1500.pdf", height = 4, width = 4)
da.plots
dev.off()

#5
pdf("Responder_Figures/BMP-like_R_vs_TSpecified_NR_da.plots.downsample1500.pdf", height = 4, width = 4)
da.plots
dev.off()



