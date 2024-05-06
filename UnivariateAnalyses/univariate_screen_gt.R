################################################
# Univariate Screening Function with and
# Without adjusting for MRD
# Global test for M5
# For T-ALL X01 Outcome Prediction
# AES
# 2023-07-10
###############################################

library(survival)
library(survminer)
library(coxphf)
library(logistf)
library(stringr)
library(writexl)
#library(janitor)
library(globaltest)

# clin.data <- annot
# M.data <- M5.AllLesions.variants2
# subsets.list <- test.list2[1:10]
# test.uni.gt <- uni_screen_gt(clin.data=annot, M.data=M5.AllLesions.variants2, subsets.list=test.list2[1:10])

uni_screen_gt <- function(clin.data, M.data, subsets.list){
  # confirm datasets are aligned
  if(!all.equal(rownames(clin.data), rownames(M.data))) stop('Clinical Data and Molecular Data are not aligned by row!')

  
  ########Unadjusted Analyses Datasests###################################
  # combine M.data with each endpoint: OS, EFS, DFS, MRDbin01
  # OS
  M.OS.df <- cbind.data.frame(clin.data$OS, clin.data$OS.status, M.data) 
  M.OS.df <- M.OS.df[which(complete.cases(M.OS.df)),]# Remove missing data (keep only complete cases)
  colnames(M.OS.df) <- c("OS", "OS.status", colnames(M.data)) # Rename the columns
  # EFS
  M.EFS.df <- cbind.data.frame(clin.data$EFS, clin.data$EFS.status, M.data) 
  M.EFS.df <- M.EFS.df[which(complete.cases(M.EFS.df)),]# Remove missing data (keep only complete cases)
  colnames(M.EFS.df) <- c("EFS", "EFS.status", colnames(M.data)) # Rename the columns
  # DFS
  M.DFS.df <- cbind.data.frame(clin.data$DFS, clin.data$DFS.status, M.data) 
  M.DFS.df <- M.DFS.df[which(complete.cases(M.DFS.df)),]# Remove missing data (keep only complete cases)
  colnames(M.DFS.df) <- c("DFS", "DFS.status", colnames(M.data)) # Rename the columns
  M.DFS.df$DFS <- ifelse(M.DFS.df$DFS==0, 0.5, M.DFS.df$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
  # MRD.bin01
  M.MRD.df <- cbind.data.frame(clin.data$MRD.bin01, M.data)
  M.MRD.df <- M.MRD.df[which(complete.cases(M.MRD.df)),] # Remove missing data (keep only complete cases)
  colnames(M.MRD.df) <- c("MRD.bin01", colnames(M.data))   # Rename the columns
  ########Adjusted Analyses Datasets#######################################
  # combine M.data with each endpoint: OS, EFS, DFS and MRD.ord to fit adjusted analysis
  # OS
  M.OS.df.adj <- cbind.data.frame(clin.data$OS, clin.data$OS.status, clin.data$MRD.ord, M.data)
  M.OS.df.adj <- M.OS.df.adj[which(complete.cases(M.OS.df.adj)),]   # Remove missing data (keep only complete cases)
  colnames(M.OS.df.adj) <- c("OS", "OS.status", "MRD.ord", colnames(M.data)) # Rename the columns
  # EFS
  M.EFS.df.adj <- cbind.data.frame(clin.data$EFS, clin.data$EFS.status, clin.data$MRD.ord, M.data)
  M.EFS.df.adj <- M.EFS.df.adj[which(complete.cases(M.EFS.df.adj)),]   # Remove missing data (keep only complete cases)
  colnames(M.EFS.df.adj) <- c("EFS", "EFS.status", "MRD.ord", colnames(M.data)) # Rename the columns
  # DFS
  M.DFS.df.adj <- cbind.data.frame(clin.data$DFS, clin.data$DFS.status, clin.data$MRD.ord, M.data)
  M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
  colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(M.data)) # Rename the columns
  M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
  
  #######Run All Analyses#######################
  message("Fitting OS Global Test")
  OS.un <- gt(Surv(M.OS.df$OS, M.OS.df$OS.status), M.OS.df[,-c(1:2)], subsets=subsets.list)
  OS.adj <- gt(Surv(OS, OS.status)~MRD.ord, ~., data=M.OS.df.adj, subsets=subsets.list)
  message("Fitting EFS Global Test")
  EFS.un <- gt(Surv(M.EFS.df$EFS, M.EFS.df$EFS.status), M.EFS.df[,-c(1:2)], subsets=subsets.list)
  EFS.adj <- gt(Surv(EFS, EFS.status)~MRD.ord, ~., data=M.EFS.df.adj, subsets=subsets.list)
  message("Fitting DFS Global Test")
  DFS.un <- gt(Surv(M.DFS.df$DFS, M.DFS.df$DFS.status), M.DFS.df[,-c(1:2)], subsets=subsets.list)
  DFS.adj <- gt(Surv(DFS, DFS.status)~MRD.ord, ~., data=M.DFS.df.adj, subsets=subsets.list)
  message("Fitting MRD Global Test")
  MRD.un <- gt(M.MRD.df$MRD.bin01, M.MRD.df[,-c(1)], subsets=subsets.list, model="logistic")
  
  
  ###### Combine Results#########################
  #if(all.equal(OS.un$MolecularFeature, OS.adj$MolecularFeature, EFS.un$MolecularFeature, EFS.adj$MolecularFeature,
  #             DFS.un$MolecularFeature, DFS.adj$MolecularFeature, MRD.un$MolecularFeature)){
  # Convert to data.frame with BH FDR
  # OS
  OS.un.df <- result(p.adjust(OS.un, method="BH"))
  colnames(OS.un.df) <- c("OS.FDR", "OS.P", "OS.Stat", "OS.Exp", "OS.Stdev", "OS.NumCov")
  pi0.hat <- min(1, 2*mean(OS.un.df$OS.P, na.rm=T))
  OS.un.df$OS.q <- pi0.hat*p.adjust(OS.un.df$OS.P, method="fdr")
  # OS.adj
  OS.adj.df <- result(p.adjust(OS.adj, method="BH"))
  colnames(OS.adj.df) <- c("OS.FDR.adj", "OS.P.adj", "OS.Stat.adj", "OS.Exp.adj", "OS.Stdev.adj", "OS.NumCov.adj")
  pi0.hat <- min(1, 2*mean(OS.adj.df$OS.P.adj, na.rm=T))
  OS.adj.df$OS.q.adj <- pi0.hat*p.adjust(OS.adj.df$OS.P.adj, method="fdr")
  
  # EFS
  EFS.un.df <- result(p.adjust(EFS.un, method="BH"))
  colnames(EFS.un.df) <- c("EFS.FDR", "EFS.P", "EFS.Stat", "EFS.Exp", "EFS.Stdev", "EFS.NumCov")
  pi0.hat <- min(1, 2*mean(EFS.un.df$EFS.P, na.rm=T))
  EFS.un.df$EFS.q <- pi0.hat*p.adjust(EFS.un.df$EFS.P, method="fdr")
  # EFS.adj
  EFS.adj.df <- result(p.adjust(EFS.adj, method="BH"))
  colnames(EFS.adj.df) <- c("EFS.FDR.adj", "EFS.P.adj", "EFS.Stat.adj", "EFS.Exp.adj", "EFS.Stdev.adj", "EFS.NumCov.adj")
  pi0.hat <- min(1, 2*mean(EFS.adj.df$EFS.P.adj, na.rm=T))
  EFS.adj.df$EFS.q.adj <- pi0.hat*p.adjust(EFS.adj.df$EFS.P.adj, method="fdr")
  
  # DFS
  DFS.un.df <- result(p.adjust(DFS.un, method="BH"))
  colnames(DFS.un.df) <- c("DFS.FDR", "DFS.P", "DFS.Stat", "DFS.Exp", "DFS.Stdev", "DFS.NumCov")
  pi0.hat <- min(1, 2*mean(DFS.un.df$DFS.P, na.rm=T))
  DFS.un.df$DFS.q <- pi0.hat*p.adjust(DFS.un.df$DFS.P, method="fdr")
  # DFS.adj
  DFS.adj.df <- result(p.adjust(DFS.adj, method="BH"))
  colnames(DFS.adj.df) <- c("DFS.FDR.adj", "DFS.P.adj", "DFS.Stat.adj", "DFS.Exp.adj", "DFS.Stdev.adj", "DFS.NumCov.adj")
  pi0.hat <- min(1, 2*mean(DFS.adj.df$DFS.P.adj, na.rm=T))
  DFS.adj.df$DFS.q.adj <- pi0.hat*p.adjust(DFS.adj.df$DFS.P.adj, method="fdr")
  
  # MRD
  MRD.un.df <- result(p.adjust(MRD.un, method="BH"))
  colnames(MRD.un.df) <- c("MRD.FDR", "MRD.P", "MRD.Stat", "MRD.Exp", "MRD.Stdev", "MRD.NumCov")
  pi0.hat <- min(1, 2*mean(MRD.un.df$MRD.P, na.rm=T))
  MRD.un.df$MRD.q <- pi0.hat*p.adjust(MRD.un.df$MRD.P, method="fdr")
  
  res.df <- cbind.data.frame(rownames(OS.un.df), OS.un.df, OS.adj.df, EFS.un.df,
                             EFS.adj.df, DFS.un.df, DFS.adj.df, MRD.un.df)
  colnames(res.df)[1] <- c("Gene")
  gt.list <- list(OS.un, OS.adj, EFS.un, EFS.adj, DFS.un, DFS.adj, MRD.un)
  names(gt.list) <- c("OS", "OS.adj", "EFS", "EFS.adj", "DFS", "DFS.adj", "MRD")
  res.list <- list(subsets.list, gt.list, res.df)
  names(res.list) <- c("Subsets", "GT.results", "Results")
  #}
  #else{stop("Data not aligned!")}
  return(res.list)
  
}

