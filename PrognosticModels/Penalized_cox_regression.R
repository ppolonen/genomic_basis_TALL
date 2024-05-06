## Script name: statistical_analysis_multivariable_model.R
##
## Purpose of script: Run glmnet, optimize alpha&lambda using CV, test using predefined data splits
##
## Author: Dr. Petri Polonen
##
## Date Created: 2023-07-26
##
## Copyright (c) Petri Polonen, 2023
## Email: petri.polonen@gmail.com
##
##*****************************************************************************************
## Running notes: Data and splits can be obtained from synapse syn54032669
##*****************************************************************************************
##
##*****************************************************************************************
##*****************************************************************************************
##*****************************************************************************************
##
#***************************** Tools and options: *****************************************
message("[STATUS] DATE RUNNING: ", Sys.Date(), appendLF = T)
ptm <- proc.time()
options(stringsAsFactors=F)
message("[STATUS] SETTING TOOLS...", appendLF = T)

# basic packages
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(require(tidyverse))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(survival))
suppressMessages(library(glmnet))

message("[STATUS] SETTING OPTIONS...", appendLF = T)
options(scipen = 5) # turns off scientific notation
options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI

#***************************** Set Directories: *******************************************
message("[STATUS] SETTING WD...", appendLF = T)

wd="/research_jude/rgs01_jude/groups/mulligrp/projects/T-ALL_GM_Teachey/common/mulligrp/petri_work/outcome/"
out="/research_jude/rgs01_jude/groups/mulligrp/projects/T-ALL_GM_Teachey/common/mulligrp/petri_work/outcome/"

setwd(wd)
#******************************************************************************************

#***************************** Input Data: ********************************************
message("[INPUT] ", appendLF = T)

load("/research_jude/rgs01_jude/groups/mulligrp/projects/T-ALL_GM_Teachey/common/biostat/anna_work/PrognosticModels/July2023FinalData/Data/Data_1309Samples.RData")
load("/research_jude/rgs01_jude/groups/mulligrp/projects/T-ALL_GM_Teachey/common/biostat/anna_work/PrognosticModels/July2023FinalData/Data/DataSplitIDs_2023-07-18.RData")

#***************************** Computations: ********************************************
message("[STATUS] WORKING...", appendLF = T)

# Clean clinical data
annot$status <- factor(annot$OS.status, levels=c(0, 1), labels=c("Alive", "Died"))
annot$MRD.ord <- as.factor(ifelse(annot$Day.29.MRD<0.1, "Negative", ifelse(annot$Day.29.MRD<10, "Low Positive", ifelse(is.na(annot$Day.29.MRD), NA, "High Positive"))))
annot$MRD.ord <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", ifelse(annot$Day.29.MRD<10, "Low Positive", "High Positive")))
annot$MRD.N <- annot$Day.29.MRD/100
annot$MRD.bin <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", "Positive"))
annot$MRD.bin01 <- ifelse(annot$MRD.bin=="Negative", 0, 1)
annot$WBC.bin <- ifelse(is.na(annot$WBC.at.Diagnosis), NA, ifelse(annot$WBC.at.Diagnosis < 200, "Low", "High"))
annot$WBC.bin01 <- ifelse(annot$WBC.bin=="Low", 0, 1)
annot$CNS.ord <- ifelse(is.na(annot$CNS.Status), NA, ifelse(annot$CNS.Status=="Unknown", NA, ifelse(annot$CNS.Status=="CNS 1", "CNS 1", ifelse(annot$CNS.Status=="CNS 2"|annot$CNS.Status=="CNS 2a"|annot$CNS.Status=="CNS 2b"|annot$CNS.Status=="CNS 2c", "CNS 2", "CNS 3"))))
annot$CNS.bin <- ifelse(is.na(annot$CNS.ord), NA, ifelse(annot$CNS.ord=="CNS 3", "3", "Other"))
annot$CNS.bin01 <- ifelse(annot$CNS.bin=="Other", 0, 1)
annot$Trt.Arm <- ifelse(annot$Treatment.Arm=="missing", NA, annot$Treatment.Arm)
annot$race <- ifelse(annot$Race=='American Indian or Alaska Native'|annot$Race=='Asian'|annot$Race=='Native Hawaiian or other Pacific Islander', "Other", annot$Race)
annot$Sex.bin01=ifelse(annot$Sex=="Male", 1, 0)
annot$DFS[annot$DFS==0]=0.5

get.concordance=function(df.features, annot, train.id.df, test.id.df, endpoint, add.clinical.features=NULL,nfold=10, cores=1, type.measure="deviance"){
  
  models=do.call(rbind, mclapply(1:100, function(i){
    
    print(i)
    
    # Step1: Train model
    ids=train.id.df[,i]
    
    if(!is.null(add.clinical.features)){
      clinical_data=annot[,add.clinical.features,drop=F]
      df.features=cbind(clinical_data, df.features)
      ids=ids[!ids%in%rownames(clinical_data)[rowSums(is.na(clinical_data))>0]] # removed NAs
    }
    
    d.train=df.features[match(ids, rownames(df.features)),,drop=F]
    
    if(sum(is.na(d.train))>0)stop("NA value in training data")
    
    # remove duplicated features/columns - same information:
    d.train=d.train[,!duplicated(as.list(d.train)),drop=F]
    
    time.train=annot[match(ids, rownames(annot)), endpoint]
    status.train=annot[match(ids, rownames(annot)), paste0(endpoint, ".status")]
    
    rm=is.na(time.train)|is.na(status.train)
    
    # add dummy
    if(dim(d.train)[2]==1){
      d.train=cbind(d.train, d.train)
    }
    
    ptm <- proc.time()
    m=suppressMessages(fun.cox.elasticnet(data = d.train[!rm,,drop=F], 
                                          time = time.train[!rm], 
                                          status = status.train[!rm],
                                          nfold=10,
                                          percentage = 0, 
                                          min.elnet = 0,
                                          max.elnet = 1,
                                          alpha_increments = 0.1,
                                          repeats = 1,
                                          cores = 1))
    
    time=signif((proc.time() - ptm)[3],1)
    message(paste0("[STATUS] ", i, " DONE, TIME: ", signif(time/60,2), "min"), appendLF = T)
    
    coef=m$coefficients
    
    # Step1: Test model
    ids2=test.id.df[,i]
    d.test=df.features[match(ids2, rownames(df.features)),,drop=F]
    time.test=annot[match(ids2, rownames(annot)), endpoint]
    status.test=annot[match(ids2, rownames(annot)), paste0(endpoint, ".status")]
    d.test.filt=d.test[,match(rownames(coef), colnames(d.test))]
    risk_patient=as.numeric(as.numeric(coef) %*% data.matrix(t(d.test.filt)))
    
    if(!all(rownames(coef)==colnames(d.test.filt)))stop("Coef not matching!")
    
    rm=is.na(time.test)|is.na(status.test)
    
    concordance <- intsurv::cIndex(time=time.test[!rm], event=status.test[!rm], risk_score=risk_patient[!rm])[1]
    data.frame("concordance"=as.numeric(concordance), "N.Coef"=length(coef), "Alpha"=m$Alpha, "Lamda"=m$Lambda)
  }, mc.cores=cores))
}

#' @param  data data.frame of variables to fit the model, NA are not permitted in glmnet
#' @param  time Survival time
#' @param  status Survival status 1 == event
#' @param  cores Number of cores to use, parallelization is done for repeats
#' @param  nfold N-Fold for crossvalidation
#' @param  min.elnet 0-1 L1+L2 regularization term minimum value (alpha in glmnet)
#' @param  max.elnet 0-1 L1+L2 regularization term maximum value (alpha in glmnet)
#' @param  max.elnet increment for alpha in glmnet
#' @param  repeats Number of repeats in glmnet
fun.cox.elasticnet=function(data, time, status, cores=8, nfold=10, min.elnet=0, max.elnet=1, alpha_increments = 0.1, repeats=100){

  # run regularized regression
  a <- seq(min.elnet, max.elnet, alpha_increments)
  
  y=Surv(time, status)
  
  formdf1 <- as.formula(paste(" ~ ", paste(colnames(data),collapse="+")))
  x=model.matrix(formdf1, data=data)
  
  # go through alpha sequence with nfolds and repeat (100 times etc) cross validation on different sets of samples
  enet.rep=do.call(rbind, parallel::mclapply(seq(repeats), function(r){
    
    s <- do.call(rbind, lapply(a, function(i){
      cv <- cv.glmnet(x=x,y=y, family = "cox", nfold = nfold, type.measure = "deviance", alpha = i, standardize=F)
      data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
    }))
    
    return(s)
  }, mc.cores=cores))
  
  # find mean alpha and lambda using all repeats
  search_m=do.call(rbind, lapply(a, function(alpha){
    cvm=mean(enet.rep$cvm[enet.rep$alpha==alpha])
    lambda=mean(enet.rep$lambda.min[enet.rep$alpha==alpha])
    data.frame(alpha, cvm, lambda)
  }))
  
  # minimum cvm error:
  cv3 <- search_m[which.min(search_m$cvm), ]
  
  print(cv3$alpha)
  print(cv3$lambda)
  print(cv3$cvm)
  
  # fit model with tuned alpha, use lambda sequence here, not single value
  md3 <- glmnet(x=x,y=y, family = "cox", alpha = cv3$alpha, standardize=F)
  
  coefficients <- coef(md3, s = cv3$lambda)[-1,,drop=F]
  coefficients=coefficients[order(abs(coefficients[,1]), decreasing = T),,drop=F]
  active_coefficients <- coefficients[,1] != 0
  
  coefficients_all=coefficients[active_coefficients,,drop=F]
  
  comprisk=data[,match(rownames(coefficients_all), colnames(data))]
  risk_patient=as.numeric(as.numeric(coefficients_all) %*% data.matrix(t(comprisk)))
  
  df=data.frame(x=risk_patient)
  ggplot(df, aes(x=x))+
    geom_density(color="lightblue", fill="lightblue")
  
  a=list(summary(coxph(y ~ PI.train, data.frame("PI.train"=risk_patient))))
  c=list(coefficients_all)
  d=list(risk_patient)
  e=list(cv3$alpha)
  f=list(cv3$lambda)
  
  out=c(a,c,d,e,f)
  names(out)=c("PI.test", "coefficients", "risk_patient", "Alpha", "Lambda")
  return(out)
  
}


M3_combinations=data.frame(M3.subtype[,!colnames(M3.subtype)%in%c("ETP-like", "TAL1 αβ-like", "TAL1 DP-like", "TLX3", "NKX2-1", "HOXA9 TCR", "KMT2A", "MLLT10", "NUP98", "NUP214")], M3.subsubtype[,colnames(M3.subsubtype)%in%c("TLX3 DP-like", "TLX3 Immature", "NKX2-1 TCR", "NKX2-1 Other")], M3.genetic.subtype, check.names = F)

#*************************** Single type **********************************

df.features=data.frame(annot[,"MRD.bin01",drop=F])
concordance.M0.bin=get.concordance(df.features, annot, train.id.df, test.id.df,add.clinical.features = "MRD.bin01", endpoint = "EFS", cores=50)

df.features=data.frame(annot[,"MRD.N",drop=F])
concordance.M0=get.concordance(df.features, annot, train.id.df, test.id.df, add.clinical.features = "MRD.N", endpoint = "EFS", cores=50)

df.features=data.frame(annot[,c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01")])
concordance.M0.allClin=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=50)

df.features=data.frame(annot[,c("MRD.bin01", "WBC.bin01", "CNS.bin01", "Sex.bin01")])
concordance.M0.allClin.bin=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.bin01", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=50)

# Clinical variables + M1
df.features=data.frame(M1.classifying.driver)
concordance.M1=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# Clinical variables + M2
df.features=data.frame(M2.ETP_status)
concordance.M2=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# Clinical variables + M3
df.features=data.frame(M3.subtype)
concordance.M3n=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

df.features=data.frame(M3.subsubtype)
concordance.M3s=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

df.features=data.frame(M3.genetic.subtype)
concordance.M3g=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# Clinical variables + M3
df.features=data.frame(M3.subtype, M3.subsubtype, M3.genetic.subtype)
concordance.M3.All=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# Clinical variables + M3 select
df.features=data.frame(M3_combinations)
concordance.M3.sel=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# Clinical variables + M4
df.features=data.frame(M4.pathway)
concordance.M4=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# Clinical variables + M5 genes
df.features=data.frame(M5.AllLesions.genes)
concordance.M5.genes=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# Clinical variables + M5 variants
df.features=data.frame(M5.AllLesions.variants)
concordance.M5.variants=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# Clinical variables + M7
df.features=data.frame(M7.IP)
concordance.M7=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

df.features=data.frame(M7.IP[,c(1,3)])
concordance.M7.wo.myeloid=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# published classifier:
df.features=data.frame(M5.AllLesions.genes[,c("KRAS", "NRAS", "FBXW7", "NOTCH1", "PTEN")])

PTEN_mutated=df.features[,c("PTEN")]>0&rowSums(df.features[,c("FBXW7", "NOTCH1")])==0
RAS_mutated=rowSums(df.features[,c("KRAS", "NRAS")])>0&rowSums(df.features[,c("FBXW7", "NOTCH1")])==0
NOTCH1=rowSums(df.features[,c("FBXW7", "NOTCH1")])>0

group=rep("Other", dim(df.features)[1])
group[PTEN_mutated]="PTEN defects"
group[RAS_mutated]="N/K-RAS mutated"
group[NOTCH1]="NOTCH1/FBXW7 mutated"

df.features=sjmisc::to_dummy(group, var.name='label')
colnames(df.features)=unique(group)
rownames(df.features)=rownames(M5.AllLesions.genes)
df.features=janitor::clean_names(df.features)

concordance.previous.group=get.concordance(df.features, annot, train.id.df, test.id.df, "EFS", add.clinical.features=c("MRD.N"),cores=100)
concordance.previous=get.concordance(df.features, annot, train.id.df, test.id.df, "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
concordance.previous.allClinical=get.concordance(df.features, annot, train.id.df, test.id.df, "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# published classifier:
df.features=data.frame(M5.AllLesions.variants[,gsub(" .*.", "", colnames(M5.AllLesions.variants))%in%c("KRAS", "NRAS", "FBXW7", "NOTCH1", "PTEN")])
concordance.previous.variants=get.concordance(df.features, annot, train.id.df, test.id.df, "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

# published classifier:
df.features=data.frame(M5.AllLesions.variants[,gsub(" .*.", "", colnames(M5.AllLesions.variants))%in%c("KRAS", "NRAS", "FBXW7", "NOTCH1", "PTEN")])
concordance.previous.variants.allClinical=get.concordance(df.features, annot, train.id.df, test.id.df, "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)

save(list=c("concordance.M0.bin", "concordance.M0", "concordance.M0.allClin", "concordance.M1", "concordance.previous", "concordance.previous.allClinical", "concordance.previous.group", "concordance.previous.variants", "concordance.previous.variants.allClinical", "concordance.M2", "concordance.M3n", "concordance.M3s", "concordance.M3g", "concordance.M3.All", "concordance.M3.sel", "concordance.M4", "concordance.M5.genes","concordance.M5.variants", "concordance.M7", "concordance.M7.wo.myeloid"), file="Concordance_M1toM7_allClinical.Rdata")

#*************************** Combinations **********************************
# bsub -P hichip -q standard -n 100 -eo multivariable.err.txt -oo multivariable.out.txt -R "rusage[mem=4000]" "Rscript run_multivariable.R"

# Clinical variables + M3 + M5
df.features=data.frame(M3.subtype,M5.AllLesions.variants)
concordance.M3n_M5=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
save(list=c("concordance.M3n_M5"), file="Concordance_combinations_allClinical_concordance.M3n_M5.Rdata")
# bsub -P hichip -q standard -n 100 -eo multivariable.err1.txt -oo multivariable.out1.txt -R "rusage[mem=2000]" "Rscript run_subtype_combinations1.R"

# Clinical variables + M3 + M5
df.features=data.frame(M3_combinations, M5.AllLesions.variants)
concordance.M3sel_M5=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
save(list=c("concordance.M3sel_M5"), file="Concordance_combinations_allClinical_concordance.M3sel_M5.Rdata")

# Clinical variables + M3sel + M1
df.features=data.frame(M3_combinations, M1.classifying.driver)
concordance.M3sel_M5=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
save(list=c("concordance.M3sel_M1"), file="Concordance_combinations_allClinical_concordance.M3sel_M1.Rdata")

# Clinical variables + M3 + M1
df.features=data.frame(M3.subtype, M1.classifying.driver)
concordance.M3sel_M5=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
save(list=c("concordance.M3sel_M1"), file="Concordance_combinations_allClinical_concordance.M3_M1.Rdata")


# Clinical variables M1 + M3 + M5
df.features=data.frame(M1.classifying.driver, M3_combinations, M5.AllLesions.variants)
concordance.M1_M3_M5=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
save(list=c("concordance.M1_M3_M5"), file="Concordance_combinations_allClinical_concordance.M1_M3_M5.Rdata")

# Clinical variables M2 + M3 + M5
df.features=data.frame(M2.ETP_status, M3_combinations, M5.AllLesions.variants)
concordance.M2_M3_M5=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
save(list=c("concordance.M2_M3_M5"), file="Concordance_combinations_allClinical_concordance.M2_M3_M5.Rdata")

# Clinical variables M4 + M3 + M5
df.features=data.frame(M4.pathway,M3_combinations, M5.AllLesions.variants)
concordance.M4_M3_M5=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
save(list=c("concordance.M4_M3_M5"), file="Concordance_combinations_allClinical_concordance.M4_M3_M5.Rdata")

# Clinical variables M1 + M2 + M3 + M4 + M5
df.features=data.frame(M1.classifying.driver,M2.ETP_status, M4.pathway, M3_combinations, M5.AllLesions.variants)
concordance.M1_M2_M3_M4_M5=get.concordance(df.features, annot, train.id.df, test.id.df, endpoint = "EFS", add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01"), cores=100)
save(list=c("concordance.M1_M2_M3_M4_M5"), file="Concordance_combinations_allClinical_concordance.M1_M2_M3_M4_M5.Rdata")

#***************************** Output Data: ***********************************************
message("[OUTPUTS] DONE", appendLF = T)

time=signif((proc.time() - ptm)[3],1)
message(paste0("[STATUS] DONE, TIME: ", time, "s"), appendLF = T)