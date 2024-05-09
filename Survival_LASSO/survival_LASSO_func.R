#### Introduction ####
# This script provides the function to run LASSO regression on the Cox
# proportional hazards model for survival analysis

library(survival)
library(glmnet)
library(selectiveInference)

# Define function to run LASSO and save results
survLASSO <- function (DEGs, featSurv.df, 
                       covariates=c(), stratSurv="ETP",
                       ETPstat=NULL, lb=NULL, ub=NULL,
                       LOOCV=TRUE, numFolds=NULL,
                       savePre="", fig.dir="./", res.dir="./", R.dir="./") {
  
  ## Create Survival object ##
  surv.obj <- Surv(time = featSurv.df$time.OS, event = featSurv.df$status.OS==1, type = "right")
  
  # Stratify survival object? (different baseline hazard)
  if (!is.null(stratSurv)) {
    surv.obj <- stratifySurv(surv.obj, featSurv.df[[stratSurv]])
  }
  
  ## Select and scale features/covariates ##
  X <- featSurv.df[,c(DEGs$gene, covariates)]
  X <- mutate(X, across(where(is.numeric),scale))
  
  # Set upper and lower bounds to inf if not provided
  if (is.null(lb) & is.null(ub)) {boundCoeff=FALSE
  } else boundCoeff=TRUE
  if (is.null(lb)) {lb <- rep(-Inf,dim(X)[2])}
  if (is.null(ub)) {ub <- rep(Inf,dim(X)[2])}
  # Determine number of folds for LOOCV
  if (LOOCV) {numFolds=dim(X)[1]}
  
  # Run LASSO using cv.glmnet cox model
  cvfit <- cv.glmnet(as.matrix(X), surv.obj, family = "cox",
                     type.measure = "deviance",
                     nfolds=numFolds, 
                     lower.limits=lb, upper.limits=ub,
                     penalty.factor=c(rep(1, ncol(X) - length(covariates)),
                                      rep(0,length(covariates))))
  
  # Get final coefficients for best fit
  beta <- coef(cvfit, x=as.matrix(X), y=surv.obj, s="lambda.min", exact=TRUE)
  # Create a df to save
  results.df <- data.frame(feature=colnames(X[,beta[,1]!=0]),coeff=beta[beta[,1]!=0,1])
  
  # Build file name based on variables  
  file.name <- paste0(savePre, "CoxLASSO")
  if (LOOCV) {
    file.name <- paste0(file.name,"_LOOCV")
  } else file.name <- paste0(file.name,"_",numFolds,"foldCV")
  if (!is.null(stratSurv)) {
    file.name <- paste0(file.name,"_",stratSurv,"strata")
  }
  if (!is.null(ETPstat)) {
    file.name <- paste0(file.name,"_",ETPstat,"only")
  }
  if (boundCoeff) {
    file.name <- paste0(file.name,"_boundCoeff")
  }
  
  # Plot CV error and save
  print(paste0(fig.dir,file.name,".png"))
  png(filename=paste0(figure.dir,file.name,".png"))
  print(plot(cvfit))
  dev.off()

  # Plot coefficient results
  #results.plot.df <- dplyr::arrange(results.df,coeff)
  results.plot.df <- dplyr::arrange(results.df[1:(nrow(results.df)-length(covar)),],coeff)
  results.plot.df$y_order <- c(1:dim(results.plot.df)[1])
  # Make Plot
  plot.lmin <- ggplot(results.plot.df, aes(x=coeff,y=reorder(feature,y_order))) +
    geom_point(size=2.0) +
    geom_col(width=0.02) +
    geom_vline(xintercept=0, linetype="dotted") +
    ylab("Feature") +
    xlab("LASSO Coefficient")
  print(paste0(fig.dir,file.name,"_lmin.pdf"))
  pdf(paste0(figure.dir,file.name,"_lmin.pdf"))
  print(plot.lmin)
  dev.off()

  # Write results to tab-delim file
  print(paste0(res.dir,file.name,"_lmin.txt"))
  write.table(results.df, file = paste0(res.dir,file.name,"_lmin.txt"), sep = "\t",
              col.names = TRUE, row.names = FALSE)

  # Save cv.glmnet object
  saveRDS(cvfit, file = paste0(R.dir,file.name,".rds"))

  # return fit
  return(list("cv.fit"=cvfit,"X"=as.matrix(X),"y"=surv.obj))
}
