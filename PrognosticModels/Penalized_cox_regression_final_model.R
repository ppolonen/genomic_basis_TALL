## Script name: Penalized_cox_regression_final_model.R
##
## Purpose of script: Fit a penalized cox regression model for full data to predict EFS
##
## Author: Dr. Petri Polonen
##
## Date Created: 2024-06-03
##
## Copyright (c) Petri Polonen, 2024
## Email: petri.polonen@gmail.com
##
##*****************************************************************************************
## Running notes: Data can be downloaded from synapse: syn54032669
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
suppressMessages(library(dplyr))
library(survival)
require(doMC)
require(glmnet)

message("[STATUS] SETTING OPTIONS...", appendLF = T)
options(scipen = 5) # turns off scientific notation
options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI

#***************************** Set Directories: *******************************************
message("[STATUS] SETTING WD...", appendLF = T)

wd=getwd()
out=""

setwd(wd)
#******************************************************************************************

#***************************** Input Data: ********************************************
message("[INPUT] ", appendLF = T)

load("MULLI-T-ALL-GM-Teachey_outcome_data_simplified_CNV_SNV.Rdata") # download syn60580346

# Clean clinical data
annot$status <- factor(annot$OS.status, levels=c(0, 1), labels=c("Alive", "Died"))
annot$MRD.ord <- as.factor(ifelse(annot$Day.29.MRD<0.1, "Negative", ifelse(annot$Day.29.MRD<10, "Low Positive", ifelse(is.na(annot$Day.29.MRD), NA, "High Positive"))))
annot$MRD.ord <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", ifelse(annot$Day.29.MRD<10, "Low Positive", "High Positive")))
annot$MRD.N <- annot$Day.29.MRD/100 # value between 0-1
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
annot=annot[lv,]

#***************************** Computations: ********************************************
message("[STATUS] WORKING...", appendLF = T)

df.features=data.frame(M3.subtype[lv,], M5.AllLesions.variants[lv,], check.names = F)

df.features=janitor::clean_names(df.features)

# Train model
add.clinical.features=c("MRD.N", "WBC.bin01", "CNS.bin01", "Sex.bin01")
ids=rownames(df.features)
endpoint="EFS"

clinical_data=annot[,add.clinical.features,drop=F]
df.features=cbind(clinical_data, df.features)
ids=ids[!ids%in%rownames(clinical_data)[rowSums(is.na(clinical_data))>0]] # removed NAs

data=df.features[match(ids, rownames(df.features)),,drop=F]

# remove duplicated columns - same information:
data=data[,!duplicated(as.list(data))]

time=annot[match(ids, rownames(annot)), endpoint]
status=annot[match(ids, rownames(annot)), paste0(endpoint, ".status")]

rm=is.na(time)|is.na(status)

formdf1 <- as.formula(paste(" ~ ", paste(colnames(data),collapse="+")))
x=model.matrix(formdf1, data=data[!rm,])
y=Surv(time[!rm], status[!rm])

registerDoMC(cores = 10)
set.seed(2023)

fit <- glmnet::cv.glmnet(x=x, y=y, family = "cox", nfolds = 10, alpha = 0.9, standardize=F, parallel = T)

# Assign lambda value manually, keeping same number of coefficients as in the train/test splits model.
coefficients <- coef(fit, s = 0.00397)[-1,,drop=F]

coefficients=coefficients[order(abs(coefficients[,1]), decreasing = T),,drop=F]
active_coefficients <- abs(coefficients[,1]) >0

coefficients_all=coefficients[active_coefficients,,drop=F]
coefficients_all=coefficients_all[order(coefficients_all, decreasing = T),,drop=F]
dim(coefficients_all)

comprisk=data[,match(rownames(coefficients_all), colnames(data))]
risk_patient=as.numeric(as.numeric(coefficients_all) %*% data.matrix(t(comprisk)))

risk_patient2=as.numeric(coefficients_all) * as.data.frame(t(t(comprisk)))
risk_patient2=rowSums(risk_patient2)

b=apply(comprisk, 1, function(v){
  v*as.numeric(coefficients_all)
})

annot.filt=annot[match(rownames(comprisk), rownames(annot)),]

concordance <- intsurv::cIndex(time=time[!rm], event=status[!rm], risk_score=risk_patient[!rm])[1]
concordance

summary(coxph(formula = y ~ risk_patient, data = data.frame(risk_patient)))

#***************************** Output Data: ***********************************************
message("[OUTPUT] ", appendLF = T)

save(fit, file = "Final_model_M3_M5.Rdata")
save(list = c("x", "y"), file = "Final_model_M3_M5_data.Rdata")

time=signif((proc.time() - ptm)[3],1)
message(paste0("[STATUS] DONE, TIME: ", time, "s"), appendLF = T)
