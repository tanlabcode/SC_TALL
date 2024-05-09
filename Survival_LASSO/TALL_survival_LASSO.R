#### Introduction ####
# This script runs penalized survival analysis using bulk RNA-seq data
# and finds optimal gene signature from DEGs determined with
# scRNA-seq analysis

#### Load Packages ####
library(tidyverse)
library(survival)
library(MASS)
library(ggplot2)
library(ggsurvfit)
library(glmnet)
source("./survival_LASSO_func.R")

# Set directory paths
figure.dir <- "./Figures/"
results.dir <- "./LASSO_results/"
R.obj.dir <- "./R_objects/"

#### Load in and Format Data ####

# BMP17
consensus.BMP.sig <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/ETP.Non.ETP.consensus.BMP.sig.mean.log2FC.0.9.rds")
consensus.Tspec.sig.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/1_RNA-seq/Figure8_NonETP_Redo/8T_OverlapBMPETP_0434/consensus.t.spec.sig.nofilter.rds")
Tspec.10 = consensus.Tspec.sig.nofilter %>% filter(avg.fc > 0.9)
# Combine BMP-spec and T-spec
consensus.DEG <- rbind(consensus.BMP.sig,Tspec.10)

# BMP surface
consensus.bmp.surface.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/BMP_surface_phenotype/consensus.bmp.surface.nofilter.rds")
consensus.bmp.surface.nofilter$gene <- consensus.bmp.surface.nofilter$coding
# Tspec Surface
consensus.t.spec.surface.nofilter <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/BMP_surface_phenotype/consensus.t.spec.surface.nofilter.rds")
consensus.t.spec.surface.nofilter$gene <- consensus.t.spec.surface.nofilter$coding
# Combine BMP-spec and T-spec
consensus.surface <- rbind(consensus.bmp.surface.nofilter[,c(-1, -3,-5)],consensus.t.spec.surface.nofilter[,c(-1, -3,-5)])
consensus.surface$gene <- consensus.surface$coding

# BMP119
BMP.119.ETP <- read_csv("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/5_Redo_survival_Analysis/Figure2/119_gene_deg.sig.csv")

# Bulk RNA-seq (0434)
bulk.object <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/Manuscript_Figures/6_Redo_survival_finalRNA/1_harmonizemtx_make_seurat/vst.no.relapse.bulk.rna.umap.20.rds")
expr.matrix <- t(as.matrix(bulk.object@assays$bulkRNA_VST@counts))

# Survival Data
patient.data <- rownames_to_column(bulk.object@meta.data,var="Patient_ID")

# create dataframe
expr.df <- rownames_to_column(as.data.frame(expr.matrix),var="Patient_ID")
# Combine to one dataframe and convert select columns to factors
full.df <- inner_join(expr.df, patient.data, by="Patient_ID") %>%
  mutate(Protocol=as.factor(Clinical.Trial.Protocol.Name),
         CNS.Status=as.factor(CNS.Status), .keep="unused")

# Convert categorical variables to numerical for LASSO
CNS.lookup <- list("CNS 1"=1, "CNS 2"=2, "CNS 2a"=2, "CNS 2b"=2, "CNS 2c"=2,
                   "CNS 3"=3, "CNS 3a"=3, "CNS 3b"=3, "CNS 3c"=3, "Unknown"=1)
full.df$CNS.Status <- sapply(full.df$CNS.Status, function(cns) return(CNS.lookup[[cns]]))
Protocol.lookup <- list("AALL03B1; AALL0434"=1, "AALL08B1; AALL0434"=2)
full.df$Protocol <- sapply(full.df$Protocol, function(p) return(Protocol.lookup[[p]]))

#### Run LASSO Function ####
#### Consensus BMP markers ####
BMP.DEG <- consensus.BMP.sig
save.pre <- "BMP17_adjAge_WBC_CNS_Protocol_"
covar <- c("dx.age", "dx.wbc", "CNS.Status", "Protocol")
# set bounds for coefficients
lower.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(0)
  } else return(-Inf) # enriched in t-spec
})
upper.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(Inf)
  } else return(0) # enriched in t-spec
})
lower.bound <- c(lower.bound,rep(-Inf,length(covar))) # no bounds for covariates
upper.bound <- c(upper.bound,rep(Inf,length(covar))) # no bounds for covariates

# Run with bounds
survLASSO.list <- survLASSO(BMP.DEG, full.df, stratSurv="ETP", ETPstat=NULL,
                            covariates=covar, 
                            lb=lower.bound, ub=upper.bound,
                            LOOCV=FALSE, numFolds=100,
                            savePre=save.pre, fig.dir=figure.dir,
                            res.dir=results.dir, R.dir=R.obj.dir)


#### Consensus Tspec markers ####
BMP.DEG <- Tspec.10
save.pre <- "Tspec10_adjAge_WBC_CNS_Protocol_"
covar <- c("dx.age", "dx.wbc", "CNS.Status", "Protocol")
# set bounds for coefficients
lower.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(0)
  } else return(-Inf) # enriched in t-spec
})
upper.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(Inf)
  } else return(0) # enriched in t-spec
})
lower.bound <- c(lower.bound,rep(-Inf,length(covar))) # no bounds for covariates
upper.bound <- c(upper.bound,rep(Inf,length(covar))) # no bounds for covariates

# Run with bounds
survLASSO.list <- survLASSO(BMP.DEG, tumor.df, stratSurv="ETP", ETPstat=NULL,
                            covariates=covar,
                            lb=lower.bound, ub=upper.bound,
                            LOOCV=FALSE, numFolds=100,
                            savePre=save.pre)


#### Consensus BMP and Tspec markers ####
BMP.DEG <- consensus.DEG
save.pre <- "BMP17andTspec10_adjAge_WBC_CNS_Protocol_"
covar <- c("dx.age", "dx.wbc", "CNS.Status", "Protocol")
# set bounds for coefficients
lower.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(0)
  } else return(-Inf) # enriched in t-spec
})
upper.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(Inf)
  } else return(0) # enriched in t-spec
})
lower.bound <- c(lower.bound,rep(-Inf,length(covar))) # no bounds for covariates
upper.bound <- c(upper.bound,rep(Inf,length(covar))) # no bounds for covariates

# Run with bounds
survLASSO.list <- survLASSO(BMP.DEG, tumor.df, stratSurv="ETP", ETPstat=NULL,
                            covariates=covar,
                            lb=lower.bound, ub=upper.bound,
                            LOOCV=FALSE, numFolds=100,
                            savePre=save.pre, fig.dir=figure.dir,
                            res.dir=results.dir, R.dir=R.obj.dir)


#### Consensus BMP surface markers ####
BMP.DEG <- consensus.bmp.surface.nofilter
save.pre <- "BMPSurface_adjAge_WBC_CNS_Protocol_"
covar <- c("dx.age", "dx.wbc", "CNS.Status", "Protocol")
# set bounds for coefficients
lower.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(0)
  } else return(-Inf) # enriched in t-spec
})
upper.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(Inf)
  } else return(0) # enriched in t-spec
})
lower.bound <- c(lower.bound,rep(-Inf,length(covar))) # no bounds for covariates
upper.bound <- c(upper.bound,rep(Inf,length(covar))) # no bounds for covariates
# Run with bounds
survLASSO.list <- survLASSO(BMP.DEG, tumor.df, stratSurv="ETP", ETPstat=NULL,
                            covariates=covar,
                            lb=lower.bound, ub=upper.bound,
                            LOOCV=FALSE, numFolds=100,
                            savePre=save.pre, fig.dir=figure.dir,
                            res.dir=results.dir, R.dir=R.obj.dir)


#### Consensus Tspec surface markers ####
BMP.DEG <- consensus.t.spec.surface.nofilter
save.pre <- "TspecSurface_adjAge_WBC_CNS_Protocol_"
covar <- c("dx.age", "dx.wbc", "CNS.Status", "Protocol")
# set bounds for coefficients
lower.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(0)
  } else return(-Inf) # enriched in t-spec
})
upper.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(Inf)
  } else return(0) # enriched in t-spec
})
lower.bound <- c(lower.bound,rep(-Inf,length(covar))) # no bounds for covariates
upper.bound <- c(upper.bound,rep(Inf,length(covar))) # no bounds for covariates

# Run with bounds
survLASSO.list <- survLASSO(BMP.DEG, tumor.df, stratSurv="ETP", ETPstat=NULL,
                            covariates=covar,
                            lb=lower.bound, ub=upper.bound,
                            LOOCV=FALSE, numFolds=100,
                            savePre=save.pre, fig.dir=figure.dir,
                            res.dir=results.dir, R.dir=R.obj.dir)


#### Consensus BMP and Tspec surface markers ####
BMP.DEG <- consensus.surface
save.pre <- "BMPandTspecSurface_adjAge_WBC_CNS_Protocol_"
covar <- c("dx.age", "dx.wbc", "CNS.Status", "Protocol")
# set bounds for coefficients
lower.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(0)
  } else return(-Inf) # enriched in t-spec
})
upper.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg.fc[i]<0) { # enriched in BMP
    return(Inf)
  } else return(0) # enriched in t-spec
})
lower.bound <- c(lower.bound,rep(-Inf,length(covar))) # no bounds for covariates
upper.bound <- c(upper.bound,rep(Inf,length(covar))) # no bounds for covariates

# Run with bounds
survLASSO.list <- survLASSO(BMP.DEG, tumor.df, stratSurv="ETP", ETPstat=NULL,
                            covariates=covar,
                            lb=lower.bound, ub=upper.bound,
                            LOOCV=FALSE, numFolds=100,
                            savePre=save.pre, fig.dir=figure.dir,
                            res.dir=results.dir, R.dir=R.obj.dir)


#### ETP only BMP markers (n=119) ####
BMP.DEG <- BMP.119.ETP
save.pre <- "BMP119_adjAge_WBC_CNS_Protocol_"
covar <- c("dx.age", "dx.wbc", "CNS.Status", "Protocol")
# set bounds for coefficients
lower.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg_log2FC[i]<0) { # enriched in BMP
    return(0)
  } else return(-Inf) # enriched in t-spec
})
upper.bound <- sapply(1:length(BMP.DEG$gene), function(i) {
  if (BMP.DEG$avg_log2FC[i]<0) { # enriched in BMP
    return(Inf)
  } else return(0) # enriched in t-spec
})
lower.bound <- c(lower.bound,rep(-Inf,length(covar))) # no bounds for covariates
upper.bound <- c(upper.bound,rep(Inf,length(covar))) # no bounds for covariates

# Run on only ETP data with LOOCV
survLASSO.list <- survLASSO(BMP.DEG, tumor.df, stratSurv=NULL, ETPstat="ETP",
                            covariates=covar,
                            lb=lower.bound, ub=upper.bound,
                            LOOCV=TRUE, numFolds=NULL,
                            #LOOCV=FALSE, numFolds=100,
                            savePre=save.pre, fig.dir=figure.dir,
                            res.dir=results.dir, R.dir=R.obj.dir)

# Run on all data with nfolds=100
survLASSO.list <- survLASSO(BMP.DEG, tumor.df, stratSurv="ETP", ETPstat=NULL,
                            covariates=covar,
                            lb=lower.bound, ub=upper.bound,
                            #LOOCV=TRUE, numFolds=NULL,
                            LOOCV=FALSE, numFolds=100,
                            savePre=save.pre, fig.dir=figure.dir,
                            res.dir=results.dir, R.dir=R.obj.dir)

