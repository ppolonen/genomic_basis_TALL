##################################################
# Univariate Firth Penalized Cox Models for 
# T-ALL X01 data
# AES
# 2023-06-30
##################################################

# load survival packages
library(survival)
library(survminer)
library(coxphf)
library(logistf)
library(stringr)
library(writexl)

# Source the univariate screening function
source("univariate_screen.R")

# Define results directory
res.dir <- getwd()

# load data
# data path on cluster
#load("/home/aseffern/PetriProject/Jan2023Data/MULLI-T-ALL-GM-Teachey_outcome_data.Rdata")
# data path locally
load("Data_1309Samples.RData")

# Create grouped variables
annot$MRD.ord <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", ifelse(annot$Day.29.MRD<10, "Low Positive", "High Positive")))
annot$MRD.bin <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", "Positive"))
annot$MRD.bin01 <- ifelse(annot$MRD.bin=="Negative", 0, 1)

### M1 Classifying Drivers
# Remove problematic characters from column names
M1.classifying.driver2 <- M1.classifying.driver
colnames(M1.classifying.driver2) <- str_replace_all(colnames(M1.classifying.driver2), "-", "_")

M1.screen.df <- uni_screen(annot, M1.classifying.driver2)
M1.screen.df$MolecularFeatureOld <- colnames(M1.classifying.driver)
M1.screen.df <- M1.screen.df[,c(1,65,2:64)]
setwd(res.dir)
write_xlsx(M1.screen.df, path="M1_Uni_Screen_2023-06-30.xlsx")

# # Factor Analysis - largest column as reference
m1.f <- as.factor(M1.classifying.driver.factor$Driver.Gene)
m1.f <- str_replace_all(m1.f, "-", "_")
levels(as.factor(m1.f))
M1.factor <- as.data.frame(m1.f)
rownames(M1.factor) <- rownames(M1.classifying.driver)
colnames(M1.factor) <- "Driver"
which.max(table(M1.factor$Driver))
M1.factor$Driver <- as.factor(M1.factor$Driver)
M1.factor$Driver <- with(M1.factor, relevel(Driver, ref="TAL1"))
M1.screen.df2 <- uni_screen(annot, M1.factor)
M1.screen.df2$MolecularFeatureOld <- str_replace_all(M1.screen.df2$MolecularFeature, "_", "-")
M1.screen.df2$MolecularFeatureOld2 <- substring(M1.screen.df2$MolecularFeatureOld, 7)
M1.screen.df2 <- M1.screen.df2[,c(1,65,66,2:64)]
setwd(res.dir)
write_xlsx(M1.screen.df2, path="M1_Uni_Screen_Factor_2023-06-30.xlsx")


### M2 ETP Status
# probably want to treat differently, combine into single factor and fit models?
M2.screen.df <- uni_screen(annot, M2.ETP_status)
setwd(res.dir)
write_xlsx(M2.screen.df, path="M2_Uni_Screen_2023-06-30.xlsx")
# Correct DFS.adj
M2.uni <- read_xlsx("M2_Uni_Screen_2023-06-30.xlsx")
# DFS Adjusted
vars.to.fix2 <- as.character(M2.uni[which(as.logical(M2.uni$DFS.WarnFlag.adj)),1])
vars.to.fix2
vars.fix2 <- M2.ETP_status[,vars.to.fix2,drop=FALSE]
M.DFS.df.adj <- cbind.data.frame(annot$DFS, annot$DFS.status, annot$MRD.ord, vars.fix2)
M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(vars.fix2)) # Rename the columns
M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
DFS.adj <- surv_analysis(vars.fix2, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=14000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("DFS", colnames(M2.uni))&grepl("adj", colnames(M2.uni)))
M2.uni2 <- M2.uni
M2.uni2[which(as.logical(M2.uni2$DFS.WarnFlag.adj)), col.ind] <- DFS.adj[,-1]
M2.uni2$DFS.FDR.adj <- p.adjust(M2.uni2$DFS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(M2.uni2$DFS.P.adj, na.rm=T))
M2.uni2$DFS.q.adj <- pi0.hat*p.adjust(M2.uni2$DFS.P.adj, method="fdr")
write.csv(M2.uni2, "M2_Uni_Screen_Corrected_2023-07-13.csv")


# Factor analysis
M2.ETP_status.factor.df <- as.data.frame(M2.ETP_status.factor)
rownames(M2.ETP_status.factor.df) <- rownames(M2.ETP_status)
colnames(M2.ETP_status.factor.df) <- "Status"
#M2.ETP_status.factor.df$Status <- as.factor(M2.ETP_status.factor.df$Status)
M2.ETP_status.factor.df$Status <- ifelse(M2.ETP_status.factor.df$Status=="Unknown"|M2.ETP_status.factor.df$Status=="missing", NA, M2.ETP_status.factor.df$Status)
M2.ETP_status.factor.df$Status <- factor(M2.ETP_status.factor.df$Status, levels=c("Non-ETP", "ETP", "Near-ETP"))
M2.screen.df2 <- uni_screen(annot, M2.ETP_status.factor.df)
setwd(res.dir)
write_xlsx(M2.screen.df2, path="M2_Uni_Screen_Factor_2023-06-30.xlsx")
# Update those that fail to converge
# DFS Adjusted
setwd(res.dir)
M2.uni <- read_xlsx("M2_Uni_Screen_Factor_2023-06-30.xlsx")
vars.to.fix2 <- M2.uni[which(as.logical(M2.uni$DFS.WarnFlag.adj)),1]
vars.to.fix2
#vars.fix2 <- M2.ETP_status.factor.df[,vars.to.fix2]
M.DFS.df.adj <- cbind.data.frame(annot$DFS, annot$DFS.status, annot$MRD.ord, M2.ETP_status.factor.df)
M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(M2.ETP_status.factor.df)) # Rename the columns
M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
DFS.adj <- surv_analysis(M2.ETP_status.factor.df, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=14000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("DFS", colnames(M2.uni))&grepl("adj", colnames(M2.uni)))
M2.uni2 <- M2.uni
M2.uni2[which(as.logical(M2.uni2$DFS.WarnFlag.adj)), col.ind] <- DFS.adj[,-1]
M2.uni2$DFS.FDR.adj <- p.adjust(M2.uni2$DFS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(M2.uni2$DFS.P.adj, na.rm=T))
M2.uni2$DFS.q.adj <- pi0.hat*p.adjust(M2.uni2$DFS.P.adj, method="fdr")
write.csv(M2.uni2, "M2_Uni_Screen_Factor_Corrected_2023-07-13.csv")


### M3 Subtype
# Convert column names to remove characters that cause problems
cn.sub <- colnames(M3.subtype)
cn.sub <- str_replace_all(cn.sub, "-", "_")
cn.sub <- str_replace_all(cn.sub, "&", "_")
cn.sub <- str_replace_all(cn.sub, "α", "a")
cn.sub <- str_replace_all(cn.sub, "β", "b")
cn.sub <- str_replace_all(cn.sub, "γ", "g")
cn.sub <- str_replace_all(cn.sub, "δ", "d")
cn.sub <- str_replace_all(cn.sub, " ", "_")
cn.sub
M3.subtype2 <- M3.subtype
colnames(M3.subtype2) <- cn.sub

M3.screen.df <- uni_screen(annot, M3.subtype2)
M3.screen.df$MolecularFeatureOld <- colnames(M3.subtype)
M3.screen.df <- M3.screen.df[,c(1,65,2:64)]
setwd(res.dir)
write_xlsx(M3.screen.df, path="M3_Uni_Screen_2023-06-30.xlsx")

# Correct 3 DFS that failed to converge
M3.uni <- read_xlsx("M3_Uni_Screen_2023-06-30.xlsx")
# DFS Adjusted
vars.to.fix2 <- unlist(M3.uni[which(as.logical(M3.uni$DFS.WarnFlag.adj)),1])
vars.to.fix2
vars.fix2 <- M3.subtype2[,vars.to.fix2,drop=FALSE]
M.DFS.df.adj <- cbind.data.frame(annot$DFS, annot$DFS.status, annot$MRD.ord, vars.fix2)
M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(vars.fix2)) # Rename the columns
M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
# set.seed(1012)
# M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, runif(1), M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
DFS.adj <- surv_analysis(vars.fix2, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=14000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("DFS", colnames(M3.uni))&grepl("adj", colnames(M3.uni)))
M3.uni2 <- M3.uni
M3.uni2[which(as.logical(M3.uni2$DFS.WarnFlag.adj)), col.ind] <- DFS.adj[,-1]
M3.uni2$DFS.FDR.adj <- p.adjust(M3.uni2$DFS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(M3.uni2$DFS.P.adj, na.rm=T))
M3.uni2$DFS.q.adj <- pi0.hat*p.adjust(M3.uni2$DFS.P.adj, method="fdr")
write.csv(M3.uni2, "M3_Uni_Screen_Corrected_2023-07-13.csv")

# Factor Analysis - largest column as reference
cn.f <- as.factor(M3.subtype.factor$Reviewed.subtype)
cn.f <- str_replace_all(cn.f, "-", "_")
cn.f <- str_replace_all(cn.f, "&", "_")
cn.f <- str_replace_all(cn.f, "α", "a")
cn.f <- str_replace_all(cn.f, "β", "b")
cn.f <- str_replace_all(cn.f, "γ", "g")
cn.f <- str_replace_all(cn.f, "δ", "d")
cn.f <- str_replace_all(cn.f, " ", "_")
levels(as.factor(cn.f))
M3.factor <- as.data.frame(cn.f)
rownames(M3.factor) <- rownames(M3.subtype)
colnames(M3.factor) <- "Subtype"
which.max(table(M3.factor$Subtype))
M3.factor$Subtype <- as.factor(M3.factor$Subtype)
M3.factor$Subtype <- with(M3.factor, relevel(Subtype, ref="TAL1_DP_like"))
M3.factor.original <- as.data.frame(as.factor(M3.subtype.factor$Reviewed.subtype))
colnames(M3.factor.original) <- c("Subtype")
M3.factor.original$Subtype <- with(M3.factor.original, relevel(Subtype, ref="TAL1 DP-like"))
M3.screen.df2 <- uni_screen(annot, M3.factor)
M3.screen.df2$MolecularFeatureOld <- levels(M3.factor.original$Subtype)[2:17]
M3.screen.df2 <- M3.screen.df2[,c(1,65,2:64)]
setwd(res.dir)
write_xlsx(M3.screen.df2, path="M3_Uni_Screen_Factor_2023-06-30.xlsx")

## M3 genetic subtype
library(janitor)
M3.genetic.subtype2 <- clean_names(M3.genetic.subtype)
rownames(M3.genetic.subtype2) <- rownames(M3.genetic.subtype)

M3g.screen.df <- uni_screen(annot, M3.genetic.subtype2)
M3g.screen.df$MolecularFeatureOld <- colnames(M3.genetic.subtype)
M3g.screen.df <- M3g.screen.df[,c(1, 65, 2:64)]
setwd(res.dir)
write_xlsx(M3g.screen.df, path="M3_genetic_Uni_Screen_2023-06-30.xlsx")

# Correct 3 DFS that failed to converge
M3.uni <- read_xlsx("M3_genetic_Uni_Screen_2023-06-30.xlsx")
# DFS Adjusted
vars.to.fix2 <- unlist(M3.uni[which(as.logical(M3.uni$DFS.WarnFlag.adj)),1])
vars.to.fix2
vars.fix2 <- M3.genetic.subtype2[,vars.to.fix2,drop=FALSE]
M.DFS.df.adj <- cbind.data.frame(annot$DFS, annot$DFS.status, annot$MRD.ord, vars.fix2)
M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(vars.fix2)) # Rename the columns
M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
# set.seed(1014)
# M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, runif(1), M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
DFS.adj <- surv_analysis(vars.fix2, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=14000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("DFS", colnames(M3.uni))&grepl("adj", colnames(M3.uni)))
M3.uni2 <- M3.uni
M3.uni2[which(as.logical(M3.uni2$DFS.WarnFlag.adj)), col.ind] <- DFS.adj[,-1]
M3.uni2$DFS.FDR.adj <- p.adjust(M3.uni2$DFS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(M3.uni2$DFS.P.adj, na.rm=T))
M3.uni2$DFS.q.adj <- pi0.hat*p.adjust(M3.uni2$DFS.P.adj, method="fdr")
write.csv(M3.uni2, "M3_genetic_Uni_Screen_Corrected_2023-07-13.csv")

## M3 subsubtype
library(janitor)
M3.subsubtype2 <- clean_names(M3.subsubtype)
rownames(M3.subsubtype2) <- rownames(M3.subsubtype)

M3s.screen.df <- uni_screen(annot, M3.subsubtype2)
M3s.screen.df$MolecularFeatureOld <- colnames(M3.subsubtype)
M3s.screen.df <- M3s.screen.df[,c(1,65,2:64)]
setwd(res.dir)
write_xlsx(M3s.screen.df, path="M3_subsubtype_Uni_Screen_2023-06-30.xlsx")

### M4
# Remove problematic characters from column names
M4.pathway2 <- clean_names(M4.pathway)
M4.data <- M4.pathway2
rownames(M4.data) <- rownames(M4.pathway)

M4.screen.df <- uni_screen(annot, M4.data)
M4.screen.df$MolecularFeatureOld <- colnames(M4.pathway)
M4.screen.df <- M4.screen.df[,c(1,65,2:64)]
setwd(res.dir)
write_xlsx(M4.screen.df, path="M4_Uni_Screen_2023-06-30.xlsx")

# Correct Adjusted DFS that failed to converge
M4.uni <- read_xlsx("M4_Uni_Screen_2023-06-30.xlsx")
# DFS Adjusted
vars.to.fix2 <- unlist(M4.uni[which(as.logical(M4.uni$DFS.WarnFlag.adj)),1])
vars.to.fix2
vars.fix2 <- M4.data[,vars.to.fix2,drop=FALSE]
M.DFS.df.adj <- cbind.data.frame(annot$DFS, annot$DFS.status, annot$MRD.ord, vars.fix2)
M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(vars.fix2)) # Rename the columns
M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
# set.seed(1018)
# M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, runif(1), M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
DFS.adj <- surv_analysis(vars.fix2, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=14000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("DFS", colnames(M4.uni))&grepl("adj", colnames(M4.uni)))
M4.uni2 <- M4.uni
M4.uni2[which(as.logical(M4.uni2$DFS.WarnFlag.adj)), col.ind] <- DFS.adj[,-1]
M4.uni2$DFS.FDR.adj <- p.adjust(M4.uni2$DFS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(M4.uni2$DFS.P.adj, na.rm=T))
M4.uni2$DFS.q.adj <- pi0.hat*p.adjust(M4.uni2$DFS.P.adj, method="fdr")
write.csv(M4.uni2, "M4_Uni_Screen_Corrected_2023-07-13.csv")

### M7 IP
M7.IP2 <- clean_names(M7.IP)
rownames(M7.IP2) <- rownames(M7.IP)

M7.screen.df <- uni_screen(annot, M7.IP2)
M7.screen.df$MolecularFeatureOld <- colnames(M7.IP)
M7.screen.df <- M7.screen.df[,c(1,65,2:64)]
setwd(res.dir)
library(writexl)
write_xlsx(M7.screen.df, path="M7_Uni_Screen_2023-06-30.xlsx")

# Correct Adjusted DFS that failed to converge
M7.uni <- read_xlsx("M7_Uni_Screen_2023-06-30.xlsx")
# DFS Adjusted
vars.to.fix2 <- unlist(M7.uni[which(as.logical(M7.uni$DFS.WarnFlag.adj)),1])
vars.to.fix2
vars.fix2 <- M7.IP2[,vars.to.fix2,drop=FALSE]
M.DFS.df.adj <- cbind.data.frame(annot$DFS, annot$DFS.status, annot$MRD.ord, vars.fix2)
M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(vars.fix2)) # Rename the columns
M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
# set.seed(1018)
# M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, runif(1), M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
DFS.adj <- surv_analysis(vars.fix2, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=14000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("DFS", colnames(M7.uni))&grepl("adj", colnames(M7.uni)))
M7.uni2 <- M7.uni
M7.uni2[which(as.logical(M7.uni2$DFS.WarnFlag.adj)), col.ind] <- DFS.adj[,-1]
M7.uni2$DFS.FDR.adj <- p.adjust(M7.uni2$DFS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(M7.uni2$DFS.P.adj, na.rm=T))
M7.uni2$DFS.q.adj <- pi0.hat*p.adjust(M7.uni2$DFS.P.adj, method="fdr")
write.csv(M7.uni2, "M7_Uni_Screen_Corrected_2023-07-13.csv")

# Fit just to factor so q is correct
# M7.f <- as.factor(M7.IP.factor$cluster.immunophenotype.lowres)
# M7.f <- str_replace_all(M7.f, "-", "_")
# levels(as.factor(M7.f))
# M7.factor <- as.data.frame(M7.f)
# rownames(M7.factor) <- rownames(M7.IP)
# colnames(M7.factor) <- "IP"
# which.max(table(M7.factor$IP))
# M7.factor$IP <- as.factor(M7.factor$IP)
# M7.screen.df2 <- uni_screen(annot, M7.factor)
# M7.screen.df2$MolecularFeatureOld <- levels(as.factor(M7.IP.factor$cluster.immunophenotype.lowres))[2:5]
# M7.screen.df2 <- M7.screen.df2[,c(1,65,2:64)]
# setwd(res.dir)
# library(writexl)
# write_xlsx(M7.screen.df2, path="M7_Uni_Screen_Factor_2023-06-30.xlsx")


#############################################
# M5 scratchwork
# Run on hpc because too large
# 
# M5.AllLesions.genes2 <- M5.AllLesions.genes
# colnames(M5.AllLesions.genes2) <- str_replace_all(colnames(M5.AllLesions.genes2), " ", "_")
# colnames(M5.AllLesions.genes2) <- str_replace_all(colnames(M5.AllLesions.genes2), ",", "_")
# 
# M5.screen.df <- uni_screen(annot, M5.AllLesions.genes2)
# M5.screen.df <- uni_screen(annot, M5.AllLesions.genes)
# 
# M5.screen.df2 <- uni_screen(annot, M5.AllLesions.variants)
# 
# M5.AllLesions.variants2 <- M5.AllLesions.variants
# colnames(M5.AllLesions.variants2) <- str_replace_all(colnames(M5.AllLesions.variants2), " ", "_")
# colnames(M5.AllLesions.variants2) <- str_replace_all(colnames(M5.AllLesions.variants2), "/", "_")
# colnames(M5.AllLesions.variants2) <- str_replace_all(colnames(M5.AllLesions.variants2), ",", "_")
# 
# M5.screen.df2 <- uni_screen(annot, M5.AllLesions.variants2)
# 
# 
# #######################################
# # gene expr univariate screen
# # Correct 12 OS.adj and 16 DFS.adj 
# # That fail to converge
# 
#source("Z:/ResearchHome/Groups/mulligrp/projects/T-ALL_GM_Teachey/common/biostat/anna_work/PrognosticModels/Code/Functions/univariate_screen_gexp_2023-06-02.R")


# OS.adj
setwd(res.dir)
gexp.uni <- read.csv("Gexp_Screen_2023-07-10.csv")
gexp.uni <- gexp.uni[,-1]

genes.to.fix <- gexp.uni[which(gexp.uni$OS.WarnFlag.adj),1]
genes.to.fix

iqr.scale.func <- function(gene.dat){
  m <- median(as.numeric(gene.dat))
  s <- IQR(gene.dat)
  scale.gene.dat <- (gene.dat-m)/s
  return(scale.gene.dat)
}
gexp.scale <- apply(gexp, 1, iqr.scale.func)
nan.vec <- colSums(is.nan(gexp.scale))
length(which(nan.vec>0)) #1029
gexp.scale.filt <- gexp.scale[,-which(nan.vec>0)]
all.equal(rownames(gexp.scale.filt), rownames(annot))

gexp.fix <- gexp.scale.filt[,genes.to.fix]

# OS Adjusted
M.OS.df.adj <- cbind.data.frame(annot$OS, annot$OS.status, annot$MRD.ord, gexp.fix)
M.OS.df.adj <- M.OS.df.adj[which(complete.cases(M.OS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.OS.df.adj) <- c("OS", "OS.status", "MRD.ord", colnames(gexp.fix)) # Rename the columns
OS.adj <- surv_analysis(gexp.fix, M.OS.df.adj, adj=TRUE, max.step=0.001, max.it=12000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("OS", colnames(gexp.uni))&grepl("adj", colnames(gexp.uni)))
gexp.uni2 <- gexp.uni
gexp.uni2[which(gexp.uni2$OS.WarnFlag.adj), col.ind] <- OS.adj[,-1]
gexp.uni2$OS.FDR.adj <- p.adjust(gexp.uni2$OS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(gexp.uni2$OS.P.adj, na.rm=T))
gexp.uni2$OS.q.adj <- pi0.hat*p.adjust(gexp.uni2$OS.P.adj, method="fdr")
#write.csv(gexp.uni2, "M5_genes_Uni_Screen_Corrected_2023-07-12.csv")

# DFS Adjusted
genes.to.fix2 <- gexp.uni[which(gexp.uni$DFS.WarnFlag.adj),1]
genes.to.fix2
gexp.fix2 <- gexp.scale.filt[,genes.to.fix2]

M.DFS.df.adj <- cbind.data.frame(annot$DFS, annot$DFS.status, annot$MRD.ord, gexp.fix2)
M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(gexp.fix2)) # Rename the columns
#M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
set.seed(16)
M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, runif(1), M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
DFS.adj <- surv_analysis(gexp.fix2, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=14000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("DFS", colnames(gexp.uni))&grepl("adj", colnames(gexp.uni)))
gexp.uni2[which(gexp.uni2$DFS.WarnFlag.adj), col.ind] <- DFS.adj[,-1]
gexp.uni2$DFS.FDR.adj <- p.adjust(gexp.uni2$DFS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(gexp.uni2$DFS.P.adj, na.rm=T))
gexp.uni2$DFS.q.adj <- pi0.hat*p.adjust(gexp.uni2$DFS.P.adj, method="fdr")
write.csv(gexp.uni2, "Gexp_Uni_Screen_Corrected_2023-07-12.csv")

