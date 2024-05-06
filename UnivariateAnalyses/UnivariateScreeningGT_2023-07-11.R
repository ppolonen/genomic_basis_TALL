##################################################
# Global Test for Gene lesions in
# T-ALL X01 data
# AES
# 2023-07-11
##################################################

# load survival packages
library(survival)
library(survminer)
library(coxphf)
library(logistf)
library(stringr)
library(writexl)
library(janitor)
library(globaltest)

# Source the univariate screening function
source("univariate_screen_gt.R")

# Define results directory
res.dir <- getwd()


# data path locally
load("Data_1309Samples.RData")

# Create grouped variables
annot$MRD.ord <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", ifelse(annot$Day.29.MRD<10, "Low Positive", "High Positive")))
annot$MRD.bin <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", "Positive"))
annot$MRD.bin01 <- ifelse(annot$MRD.bin=="Negative", 0, 1)

# M5 Cleaning
M5.AllLesions.genes2 <- clean_names(M5.AllLesions.genes)
rownames(M5.AllLesions.genes2) <- rownames(M5.AllLesions.genes)
M5.AllLesions.variants2 <- clean_names(M5.AllLesions.variants)
rownames(M5.AllLesions.variants2) <- rownames(M5.AllLesions.variants)

var.names <- cbind.data.frame(colnames(M5.AllLesions.variants), colnames(M5.AllLesions.variants2))
colnames(var.names) <- c("Original", "Clean")

gene.names <- cbind.data.frame(colnames(M5.AllLesions.genes), colnames(M5.AllLesions.genes2))
colnames(gene.names) <- c("Original", "Clean")

# Make gene/variant list
variant <- colnames(M5.AllLesions.variants2)
gene2 <- str_remove_all(string=colnames(M5.AllLesions.variants), pattern=" .*")
gene <- make_clean_names(gene2, allow_dupes = TRUE)
gene.var.df <- cbind.data.frame(gene, variant)
unique.gene <- unique(gene)
length(which(unique.gene %in% colnames(M5.AllLesions.genes2)))/length(unique.gene)
unique.gene[-which(unique.gene %in% colnames(M5.AllLesions.genes2))]
# Is in genes, just named slightly different
# replace
gene.var.df[grep("cn_vgains", gene.var.df$gene),1] <- "cn_vgains_8_10_19"
colnames(M5.AllLesions.genes2)[-which(colnames(M5.AllLesions.genes2) %in% unique.gene)]
colnames(M5.AllLesions.genes)[-which(colnames(M5.AllLesions.genes2) %in% unique.gene)]

# Make list
gene.list <- split(gene.var.df$variant, gene.var.df$gene)

# Analyze
test.uni.gt <- uni_screen_gt(clin.data=annot, M.data=M5.AllLesions.variants2, subsets.list=gene.list[1:10])
M5.uni.gt <- uni_screen_gt(clin.data=annot, M.data=M5.AllLesions.variants2, subsets.list=gene.list)
# Map to old M5 gene names
res <- M5.uni.gt$Results
gene.names2 <- gene.names
colnames(gene.names2)[2] <- "Gene"
res2 <- dplyr::full_join(res, gene.names2, by="Gene")
res3 <- res2[,c(1,51,2:50)]

# Save results
save.image("M5_GlobalTest.RData")
#load("Z:/ResearchHome/Groups/mulligrp/projects/T-ALL_GM_Teachey/common/biostat/anna_work/PrognosticModels/M5_GlobalTest.RData")
write.csv(res3, file=paste0(res.dir, "M5_gt_Uni_Screen_2023-07-11.csv"))
save(M5.uni.gt, gene.list, var.names, gene.names, res3, file=paste0(res.dir, "M5Screen.RData"))

# Review Results and run univariate screen on variants from significant genes
length(which(res3$OS.q<0.1))
os.sig.genes <- res3[which(res3$OS.q<0.1),]$Gene
os.adj.sig.genes <- res3[which(res3$OS.q.adj<0.1),]$Gene
efs.sig.genes <- res3[which(res3$EFS.q<0.1),]$Gene
efs.adj.sig.genes <- res3[which(res3$EFS.q.adj<0.1),]$Gene
dfs.sig.genes <- res3[which(res3$DFS.q<0.1),]$Gene
dfs.adj.sig.genes <- res3[which(res3$DFS.q.adj<0.1),]$Gene
mrd.sig.genes <- res3[which(res3$MRD.q<0.1),]$Gene
sig.genes.comb <- c(os.sig.genes, os.adj.sig.genes, efs.sig.genes, efs.adj.sig.genes, dfs.sig.genes,
                    dfs.adj.sig.genes, mrd.sig.genes)
sig.genes.un <- unique(sig.genes.comb)
sig.gene.list <- gene.list[names(gene.list) %in% sig.genes.un]
sig.variants <- unlist(sig.gene.list)
# Run univariate screen on filtered variant list
M5.vars.sig <- M5.AllLesions.variants2[,which(colnames(M5.AllLesions.variants2) %in% sig.variants)]
M5.vars.sig.oldnames <- colnames(M5.AllLesions.variants[,which(colnames(M5.AllLesions.variants2) %in% sig.variants)])
source("univariate_screen.R")
M5.screen.df <- uni_screen(annot, M5.vars.sig)
M5.screen.df$MolecularFeatureOld <- M5.vars.sig.oldnames
M5.screen.df <- M5.screen.df[,c(1,65,2:64)]
setwd(res.dir)
write_xlsx(M5.screen.df, path="M5v_Uni_Screen_2023-07-12.xlsx")

# Create final excel results of sig genes + variant analysis
M5v.OS <- M5.screen.df[unique(grep(paste(os.sig.genes, collapse="|"), M5.screen.df$MolecularFeature)), c(1,2,grep("OS", colnames(M5.screen.df)))]
M5v.OS <- M5v.OS[,-grep(paste(c("OS.FDR", "OS.q", "OS.FDR.adj", "OS.q.adj"), collapse="|"), colnames(M5v.OS))]
M5v.EFS <- M5.screen.df[unique(grep(paste(efs.sig.genes, collapse="|"), M5.screen.df$MolecularFeature)), c(1,2,grep("EFS", colnames(M5.screen.df)))]
M5v.EFS <- M5v.EFS[,-grep(paste(c("EFS.FDR", "EFS.q", "EFS.FDR.adj", "EFS.q.adj"), collapse="|"), colnames(M5v.EFS))]
M5v.DFS <- M5.screen.df[unique(grep(paste(dfs.sig.genes, collapse="|"), M5.screen.df$MolecularFeature)), c(1,2,grep("DFS", colnames(M5.screen.df)))]
M5v.DFS <- M5v.DFS[,-grep(paste(c("DFS.FDR", "DFS.q", "DFS.FDR.adj", "DFS.q.adj"), collapse="|"), colnames(M5v.DFS))]
M5v.MRD <- M5.screen.df[unique(grep(paste(mrd.sig.genes, collapse="|"), M5.screen.df$MolecularFeature)), c(1,2,grep("MRD", colnames(M5.screen.df)))]
M5v.MRD <- M5v.MRD[,-grep(paste(c("MRD.FDR", "MRD.q", "MRD.FDR.adj", "MRD.q.adj"), collapse="|"), colnames(M5v.MRD))]
tab.list <- list(res3, M5v.OS, M5v.EFS, M5v.DFS, M5v.MRD)
names(tab.list) <- c("GlobalTest", "OS", "EFS", "DFS", "MRD")
setwd(res.dir)
write_xlsx(tab.list, "M5GlobalTestAndVariantScreen_2023-07-12.xlsx")

###################################################
# Update with P-value threshold instead of q-value
# 2023-07-13
###################################################
# Review Results and run univariate screen on variants from significant genes
length(which(res3$OS.P<0.1))
os.sig.genes <- res3[which(res3$OS.P<0.1),]$Gene
os.adj.sig.genes <- res3[which(res3$OS.P.adj<0.1),]$Gene
efs.sig.genes <- res3[which(res3$EFS.P<0.1),]$Gene
efs.adj.sig.genes <- res3[which(res3$EFS.P.adj<0.1),]$Gene
dfs.sig.genes <- res3[which(res3$DFS.P<0.1),]$Gene
dfs.adj.sig.genes <- res3[which(res3$DFS.P.adj<0.1),]$Gene
mrd.sig.genes <- res3[which(res3$MRD.P<0.1),]$Gene
sig.genes.comb <- c(os.sig.genes, os.adj.sig.genes, efs.sig.genes, efs.adj.sig.genes, dfs.sig.genes,
                    dfs.adj.sig.genes, mrd.sig.genes)
sig.genes.un <- unique(sig.genes.comb)
sig.gene.list <- gene.list[names(gene.list) %in% sig.genes.un]
sig.variants <- unlist(sig.gene.list)
# Run univariate screen on filtered variant list
M5.vars.sig <- M5.AllLesions.variants2[,which(colnames(M5.AllLesions.variants2) %in% sig.variants)]
M5.vars.sig.oldnames <- colnames(M5.AllLesions.variants[,which(colnames(M5.AllLesions.variants2) %in% sig.variants)])
source("Z:/ResearchHome/Groups/mulligrp/projects/T-ALL_GM_Teachey/common/biostat/anna_work/PrognosticModels/Code/Functions/univariate_screen.R")
M5.screen.df <- uni_screen(annot, M5.vars.sig)
M5.screen.df$MolecularFeatureOld <- M5.vars.sig.oldnames
M5.screen.df <- M5.screen.df[,c(1,65,2:64)]
M5.uni <- M5.screen.df
v1 <- c("Model", paste0("Number Fail to Converge out of ", nrow(M5.uni)))
v2 <- c("OS", length(M5.uni$MolecularFeature[which(M5.uni$OS.WarnFlag=="TRUE")]))
v3 <- c("OS.adj", length(M5.uni$MolecularFeature[which(M5.uni$OS.WarnFlag.adj=="TRUE")]))
v4 <- c("EFS", length(M5.uni$MolecularFeature[which(M5.uni$EFS.WarnFlag=="TRUE")]))
v5 <- c("EFS.adj", length(M5.uni$MolecularFeature[which(M5.uni$EFS.WarnFlag.adj=="TRUE")]))
v6 <- c("DFS", length(M5.uni$MolecularFeature[which(M5.uni$DFS.WarnFlag=="TRUE")]))
v7<- c("DFS.adj", length(M5.uni$MolecularFeature[which(M5.uni$DFS.WarnFlag.adj=="TRUE")]))
v8 <- c("MRD", length(M5.uni$MolecularFeature[which(M5.uni$MRD.WarnFlag=="TRUE")]))
df <- rbind.data.frame(v2, v3, v4, v5, v6, v7, v8)
colnames(df) <- v1
df
# Correct DFS for some
vars.to.fix2 <- unlist(M5.uni[which(as.logical(M5.uni$DFS.WarnFlag.adj)),1])
vars.to.fix2
vars.fix2 <- M5.vars.sig[,vars.to.fix2,drop=FALSE]
M.DFS.df.adj <- cbind.data.frame(annot$DFS, annot$DFS.status, annot$MRD.ord, vars.fix2)
M.DFS.df.adj <- M.DFS.df.adj[which(complete.cases(M.DFS.df.adj)),]   # Remove missing data (keep only complete cases)
colnames(M.DFS.df.adj) <- c("DFS", "DFS.status", "MRD.ord", colnames(vars.fix2)) # Rename the columns
M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, 0.5, M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
# set.seed(1018)
# M.DFS.df.adj$DFS <- ifelse(M.DFS.df.adj$DFS==0, runif(1), M.DFS.df.adj$DFS) # Adjust the 0 survival times to 0.5 so models can be fit
DFS.adj <- surv_analysis(vars.fix2, M.DFS.df.adj, adj=TRUE, max.step=0.001, max.it=14000)
# Add estimates and P-value to spreadsheet, then recalculate FDR and q because based on all p-values
col.ind <- which(grepl("DFS", colnames(M5.uni))&grepl("adj", colnames(M5.uni)))
M5.uni2 <- M5.uni
M5.uni2[which(as.logical(M5.uni2$DFS.WarnFlag.adj)), col.ind] <- DFS.adj[,-1]
M5.uni2$DFS.FDR.adj <- p.adjust(M5.uni2$DFS.P.adj, method="fdr")
pi0.hat <- min(1, 2*mean(M5.uni2$DFS.P.adj, na.rm=T))
M5.uni2$DFS.q.adj <- pi0.hat*p.adjust(M5.uni2$DFS.P.adj, method="fdr")
setwd(res.dir)
write_xlsx(M5.uni2, path="M5v_Uni_Screen_PThreshold_2023-07-13.xlsx")

# Create final excel results of sig genes + variant analysis
M5.screen.df <- M5.uni2
M5v.OS <- M5.screen.df[unique(grep(paste(os.sig.genes, collapse="|"), M5.screen.df$MolecularFeature)), c(1,2,grep("OS", colnames(M5.screen.df)))]
M5v.OS <- M5v.OS[,-grep(paste(c("OS.FDR", "OS.q", "OS.FDR.adj", "OS.q.adj"), collapse="|"), colnames(M5v.OS))]
M5v.EFS <- M5.screen.df[unique(grep(paste(efs.sig.genes, collapse="|"), M5.screen.df$MolecularFeature)), c(1,2,grep("EFS", colnames(M5.screen.df)))]
M5v.EFS <- M5v.EFS[,-grep(paste(c("EFS.FDR", "EFS.q", "EFS.FDR.adj", "EFS.q.adj"), collapse="|"), colnames(M5v.EFS))]
M5v.DFS <- M5.screen.df[unique(grep(paste(dfs.sig.genes, collapse="|"), M5.screen.df$MolecularFeature)), c(1,2,grep("DFS", colnames(M5.screen.df)))]
M5v.DFS <- M5v.DFS[,-grep(paste(c("DFS.FDR", "DFS.q", "DFS.FDR.adj", "DFS.q.adj"), collapse="|"), colnames(M5v.DFS))]
M5v.MRD <- M5.screen.df[unique(grep(paste(mrd.sig.genes, collapse="|"), M5.screen.df$MolecularFeature)), c(1,2,grep("MRD", colnames(M5.screen.df)))]
M5v.MRD <- M5v.MRD[,-grep(paste(c("MRD.FDR", "MRD.q", "MRD.FDR.adj", "MRD.q.adj"), collapse="|"), colnames(M5v.MRD))]
tab.list <- list(res3, M5v.OS, M5v.EFS, M5v.DFS, M5v.MRD)
names(tab.list) <- c("GlobalTest", "OS", "EFS", "DFS", "MRD")
setwd(res.dir)
write_xlsx(tab.list, "M5GlobalTestAndVariantScreenP_2023-07-13.xlsx")

