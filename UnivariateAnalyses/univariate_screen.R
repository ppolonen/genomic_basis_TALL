################################################
# Univariate Screening Function with and
# Without adjusting for MRD
# For T-ALL X01 Outcome Prediction
# AES
# 2023-01-26
###############################################

library(survival)
library(survminer)
library(coxphf)
library(logistf)

#clin.data <- annot
#M.data <- M3.factor.ss

uni_screen <- function(clin.data, M.data){
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
  message("Fitting OS Cox Models")
  OS.un <- surv_analysis(M.data, M.OS.df, adj=FALSE, max.step=0.01, max.it=500)
  OS.adj <- surv_analysis(M.data, M.OS.df.adj, adj=TRUE, max.step=0.01, max.it=500)
  message("Fitting EFS Cox Models")
  EFS.un <- surv_analysis(M.data, M.EFS.df, adj=FALSE, max.step=0.01, max.it=500)
  EFS.adj <- surv_analysis(M.data, M.EFS.df.adj, adj=TRUE, max.step=0.01, max.it=500)
  message("Fitting DFS Cox Models")
  DFS.un <- surv_analysis(M.data, M.DFS.df, adj=FALSE, max.step=0.01, max.it=500)
  DFS.adj <- surv_analysis(M.data, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=1000)
  message("Fitting MRD Logistic Regression")
  MRD.un <- lr_analysis(M.data, M.MRD.df)
  
  
  ###### Combine Results#########################
  #if(all.equal(OS.un$MolecularFeature, OS.adj$MolecularFeature, EFS.un$MolecularFeature, EFS.adj$MolecularFeature,
  #             DFS.un$MolecularFeature, DFS.adj$MolecularFeature, MRD.un$MolecularFeature)){
    res.df <- cbind.data.frame(OS.un, OS.adj[,-1], EFS.un[,-1], EFS.adj[,-1], DFS.un[,-1], DFS.adj[,-1], MRD.un[,-1])
  #}
  #else{stop("Data not aligned!")}
  return(res.df)
  
}


#######################################################################
# Firth-Penalized Cox Model Screening with and without MRD adjustment
######################################################################
#M.data <- M1.classifying.driver
#M.surv.data <- M.OS.df.adj
surv_analysis <- function(M.data, M.surv.data, adj=FALSE, max.step=0.1, max.it=100){
  if(adj){
    M1.OS.list <- list()
    M1.OS.res.list <- list()
    for(i in 1:ncol(M.data)){
      temp.driver <- colnames(M.data)[i] # Extract a classifying driver
      temp.data <- M.surv.data[,c(1,2, 3, grep(paste("^", temp.driver, "$", sep=""), colnames(M.surv.data)))] # build a temporary dataset
      colnames(temp.data)[1:2] <- c("time", "status")
      test.cox <- tryCatch(coxphf(Surv(time, status)~., data=temp.data, maxstep=max.step, maxit=max.it), warning=function(w) w)
      hasWarn <- is(test.cox, "warning")
      temp.cox <- coxphf(Surv(time, status)~., data=temp.data, maxstep=max.step, maxit=max.it)
      M1.OS.list[[i]] <- temp.cox
      if(length(temp.cox$coefficients)==3){
        M1.OS.res.list[[i]]<- unname(c(temp.driver, temp.cox$coefficients[3], temp.cox$ci.lower[3], temp.cox$ci.upper[3], qchisq(1-temp.cox$prob[3],1), temp.cox$prob[3], hasWarn))
      }
      else{
        temp.df <- data.frame()
        for(j in 1:length(temp.cox$coefficients)){
          v1 <- unname(c(names(temp.cox$coefficients)[j], temp.cox$coefficients[j], temp.cox$ci.lower[j], temp.cox$ci.upper[j], qchisq(1-temp.cox$prob[j],1), temp.cox$prob[j], hasWarn))
          temp.df <- rbind.data.frame(temp.df, v1)
        }
        M1.OS.res.list[[i]] <- temp.df
      }
    }
    M1.OS.res.df <- do.call(rbind.data.frame, M1.OS.res.list)
    colnames(M1.OS.res.df) <- c("MolecularFeature", "Coef", "Lower95", "Upper95", "Stat", "P", "WarnFlag")
    if(any(grepl("MRD", M1.OS.res.df$MolecularFeature))){ 
      M1.OS.res.df <- M1.OS.res.df[-grep("MRD", M1.OS.res.df$MolecularFeature),]}
    else{}
    M1.OS.res.df$Coef <- as.numeric(M1.OS.res.df$Coef)
    M1.OS.res.df$Lower95 <- as.numeric(M1.OS.res.df$Lower95)
    M1.OS.res.df$Upper95 <- as.numeric(M1.OS.res.df$Upper95)
    M1.OS.res.df$Stat <- as.numeric(M1.OS.res.df$Stat)
    M1.OS.res.df$P <- as.numeric(M1.OS.res.df$P)
    M1.OS.res.df$fdr <- p.adjust(M1.OS.res.df$P)
    M1.OS.res.df$HR <- exp(as.numeric(M1.OS.res.df$Coef))
    # Find q-values
    pi0.hat <- min(1, 2*mean(M1.OS.res.df$P, na.rm=T))
    M1.OS.res.df$q <- pi0.hat*p.adjust(M1.OS.res.df$P, method="fdr")
    col.order <- c("MolecularFeature", "Coef", "HR", "Lower95", "Upper95", "Stat", "P", "fdr", "q", "WarnFlag")
    M1.OS.res.df <- M1.OS.res.df[,col.order]
    #View(M1.OS.res.df)
    colnames(M1.OS.res.df) <- c("MolecularFeature", paste(colnames(M.surv.data)[1],"Coef", "adj", sep="."),
                                paste(colnames(M.surv.data)[1],"HR", "adj", sep="."),
                                paste(colnames(M.surv.data)[1],"Lower95", "adj", sep="."), 
                                paste(colnames(M.surv.data)[1],"Upper95", "adj", sep="."),
                                paste(colnames(M.surv.data)[1],"Stat", "adj", sep="."),
                                paste(colnames(M.surv.data)[1],"P", "adj", sep="."),
                                paste(colnames(M.surv.data)[1],"FDR", "adj", sep="."),
                                paste(colnames(M.surv.data)[1],"q", "adj", sep="."),
                                paste(colnames(M.surv.data)[1],"WarnFlag", "adj", sep="."))
  }
  else{
    M1.OS.list <- list()
    M1.OS.res.list <- list()
    for(i in 1:ncol(M.data)){
      temp.driver <- colnames(M.data)[i] # Extract a molecular feature
      temp.data <- M.surv.data[,c(1,2,grep(paste("^", temp.driver, "$", sep=""), colnames(M.surv.data)))] # build a temporary dataset
      colnames(temp.data)[1:2] <- c("time", "status")
      test.cox <- tryCatch(coxphf(Surv(time, status)~., data=temp.data, maxstep=max.step, maxit=max.it), warning=function(w) w)
      hasWarn <- is(test.cox, "warning")
      temp.cox <- coxphf(Surv(time, status)~., data=temp.data, maxstep=max.step, maxit=max.it)
      M1.OS.list[[i]] <- temp.cox
      if(length(temp.cox$coefficients)==1){
        M1.OS.res.list[[i]]<- c(temp.driver, temp.cox$coefficients, temp.cox$ci.lower, temp.cox$ci.upper, temp.cox$coefficients/sqrt(temp.cox$var), temp.cox$prob, hasWarn)
      }
      else{
        temp.df <- data.frame()
        for(j in 1:length(temp.cox$coefficients)){
          v1 <- c(names(temp.cox$coefficients)[j], temp.cox$coefficients[j], temp.cox$ci.lower[j], temp.cox$ci.upper[j], temp.cox$coefficients[j]/sqrt(temp.cox$var[j,j]), temp.cox$prob[j], hasWarn)
          temp.df <- rbind.data.frame(temp.df, v1)
          }
        M1.OS.res.list[[i]] <- temp.df
      }
    }
    M1.OS.res.df <- do.call(rbind.data.frame, M1.OS.res.list)
    colnames(M1.OS.res.df) <- c("MolecularFeature", "Coef", "Lower95", "Upper95", "Stat", "P", "WarnFlag")
    M1.OS.res.df$Coef <- as.numeric(M1.OS.res.df$Coef)
    M1.OS.res.df$Lower95 <- as.numeric(M1.OS.res.df$Lower95)
    M1.OS.res.df$Upper95 <- as.numeric(M1.OS.res.df$Upper95)
    M1.OS.res.df$Stat <- as.numeric(M1.OS.res.df$Stat)
    M1.OS.res.df$P <- as.numeric(M1.OS.res.df$P)
    M1.OS.res.df$fdr <- p.adjust(M1.OS.res.df$P)
    M1.OS.res.df$HR <- exp(as.numeric(M1.OS.res.df$Coef))
    # Find q-values
    pi0.hat <- min(1, 2*mean(M1.OS.res.df$P, na.rm=T))
    M1.OS.res.df$q <- pi0.hat*p.adjust(M1.OS.res.df$P, method="fdr")
    col.order <- c("MolecularFeature", "Coef", "HR", "Lower95", "Upper95", "Stat", "P", "fdr", "q", "WarnFlag")
    M1.OS.res.df <- M1.OS.res.df[,col.order]
    #View(M1.OS.res.df)
    colnames(M1.OS.res.df) <- c("MolecularFeature", paste(colnames(M.surv.data)[1],"Coef", sep="."),
                                paste(colnames(M.surv.data)[1],"HR", sep="."),
                                paste(colnames(M.surv.data)[1],"Lower95", sep="."), 
                                paste(colnames(M.surv.data)[1],"Upper95", sep="."),
                                paste(colnames(M.surv.data)[1],"Stat", sep="."),
                                paste(colnames(M.surv.data)[1],"P", sep="."),
                                paste(colnames(M.surv.data)[1],"FDR", sep="."),
                                paste(colnames(M.surv.data)[1],"q", sep="."),
                                paste(colnames(M.surv.data)[1],"WarnFlag", sep="."))
  }
  
  return(M1.OS.res.df)
}

###########################################################
# Firth-penalized Logistic Regression for MRD
###########################################################

lr_analysis <- function(M.data, M.MRD.data){
  M1.MRD.list <- list()
  M1.MRD.res.list <- list()
  for(i in 1:ncol(M.data)){
    temp.driver <- colnames(M.data)[i] # Extract a molecular feature
    temp.data <- M.MRD.data[,c(1,grep(paste("^", temp.driver, "$", sep=""), colnames(M.MRD.data)))] # build a temporary dataset
    test.lr <- tryCatch(logistf(MRD.bin01~., data=temp.data, maxstep=max.step, maxit=max.it), warning=function(w) w)
    hasWarn <- is(test.lr, "warning")
    temp.lr <- logistf(MRD.bin01~., data=temp.data, maxstep=max.step, maxit=max.it)
    M1.MRD.list[[i]] <- temp.lr
    if(length(temp.lr$coefficients)==2){
      M1.MRD.res.list[[i]]<- unname(c(temp.driver, temp.lr$coefficients[2], temp.lr$ci.lower[2], temp.lr$ci.upper[2], 2*(temp.lr$loglik[1]-temp.lr$loglik[2]), temp.lr$prob[2], hasWarn))
    }
    else{
      temp.df <- data.frame()
      for(j in 2:length(temp.lr$coefficients)){
        v1 <- unname(c(names(temp.lr$coefficients)[j], temp.lr$coefficients[j], temp.lr$ci.lower[j], temp.lr$ci.upper[j], 2*(temp.lr$loglik[1]-temp.lr$loglik[2]), temp.lr$prob[j], hasWarn))
        temp.df <- rbind.data.frame(temp.df, v1)
      }
      M1.MRD.res.list[[i]] <- temp.df
    }
  }
  M1.MRD.res.df <- do.call(rbind.data.frame, M1.MRD.res.list)
  colnames(M1.MRD.res.df) <- c("MolecularFeature", "Coef", "Lower95", "Upper95", "Stat", "P", "WarnFlag")
  M1.MRD.res.df$Coef <- as.numeric(M1.MRD.res.df$Coef)
  M1.MRD.res.df$Lower95 <- as.numeric(M1.MRD.res.df$Lower95)
  M1.MRD.res.df$Upper95 <- as.numeric(M1.MRD.res.df$Upper95)
  M1.MRD.res.df$Stat <- as.numeric(M1.MRD.res.df$Stat)
  M1.MRD.res.df$P <- as.numeric(M1.MRD.res.df$P)
  M1.MRD.res.df$fdr <- p.adjust(M1.MRD.res.df$P)
  M1.MRD.res.df$OR <- exp(as.numeric(M1.MRD.res.df$Coef))
  # Find q-values
  pi0.hat <- min(1, 2*mean(M1.MRD.res.df$P, na.rm=T))
  M1.MRD.res.df$q <- pi0.hat*p.adjust(M1.MRD.res.df$P, method="fdr")
  col.order <- c("MolecularFeature", "Coef", "OR", "Lower95", "Upper95", "Stat", "P", "fdr", "q", "WarnFlag")
  M1.MRD.res.df <- M1.MRD.res.df[,col.order]
  #View(M1.MRD.res.df)
  colnames(M1.MRD.res.df) <- c("MolecularFeature", "MRD.Coef", "MRD.OR", "MRD.Lower95",
                               "MRD.Upper95", "MRD.Stat", "MRD.P", "MRD.FDR",
                               "MRD.q", "MRD.WarnFlag")
  return(M1.MRD.res.df)
}
