##################################################
# Univariate Firth Penalized Cox Models for 
# T-ALL X01 data
# Pairwise MRD comparisons (Negative, RD, IF)
# AES
# 2023-07-16
##################################################

# load survival packages
library(survival)
library(survminer)
library(coxphf)
library(logistf)
library(stringr)
library(writexl)

# Source the univariate screening function
source("univariate_screen_pairwiseMRD.R")
# Source global test univariate screening function
source("univariate_screen_pairwiseMRD_gt.R")


# Define results directory
res.dir <- getwd()

# load data
load("\Data_1309Samples.RData")

# Create grouped variables
annot$MRD.ord <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", ifelse(annot$Day.29.MRD<10, "Low Positive", "High Positive")))
annot$MRD.bin <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", "Positive"))
annot$MRD.bin01 <- ifelse(annot$MRD.bin=="Negative", 0, 1)
# Define MRD Negative, Relapse Disease (RD), and induction failure (IF)
annot$MRD.fac <- ifelse(is.na(annot$Day.29.MRD), NA, ifelse(annot$Day.29.MRD<0.1, "Negative", ifelse(annot$Day.29.MRD<10, "RD", "IF")))
annot$MRD.NegRD <- ifelse(annot$MRD.fac=="RD", 1, ifelse(annot$MRD.fac=="Negative", 0, NA))
annot$MRD.NegIF <- ifelse(annot$MRD.fac=="IF", 1, ifelse(annot$MRD.fac=="Negative", 0, NA))
annot$MRD.RDIF <- ifelse(annot$MRD.fac=="RD", 0, ifelse(annot$MRD.fac=="IF", 1, NA))

#######M1 Classifying Drivers################
# Remove problematic characters from column names
M1.classifying.driver2 <- M1.classifying.driver
colnames(M1.classifying.driver2) <- str_replace_all(colnames(M1.classifying.driver2), "-", "_")

M1.screen.df <- uni_screen(annot, M1.classifying.driver2)
M1.screen.df$MolecularFeatureOld <- colnames(M1.classifying.driver)
M1.screen.df <- M1.screen.df[,c(1,29,2:28)]
setwd(res.dir)
write_xlsx(M1.screen.df, path="M1_PairwiseMRD_Screen_2023-07-16.xlsx")

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
M1.screen.df2 <- M1.screen.df2[,c(1,29, 30, 2:28)]
setwd(res.dir)
write_xlsx(M1.screen.df2, path="M1_PairwiseMRD_Screen_Factor_2023-07-16.xlsx")

### M2 ETP Status
# probably want to treat differently, combine into single factor and fit models?
M2.screen.df <- uni_screen(annot, M2.ETP_status)
setwd(res.dir)
write_xlsx(M2.screen.df, path="M2_PairwiseMRD_Screen_2023-07-16.xlsx")

# Factor analysis
M2.ETP_status.factor.df <- as.data.frame(M2.ETP_status.factor)
rownames(M2.ETP_status.factor.df) <- rownames(M2.ETP_status)
colnames(M2.ETP_status.factor.df) <- "Status"
#M2.ETP_status.factor.df$Status <- as.factor(M2.ETP_status.factor.df$Status)
M2.ETP_status.factor.df$Status <- ifelse(M2.ETP_status.factor.df$Status=="Unknown"|M2.ETP_status.factor.df$Status=="missing", NA, M2.ETP_status.factor.df$Status)
M2.ETP_status.factor.df$Status <- factor(M2.ETP_status.factor.df$Status, levels=c("Non-ETP", "ETP", "Near-ETP"))
M2.screen.df2 <- uni_screen(annot, M2.ETP_status.factor.df)
setwd(res.dir)
write_xlsx(M2.screen.df2, path="M2_PairwiseMRD_Screen_Factor_2023-07-16.xlsx")

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
M3.screen.df <- M3.screen.df[,c(1,29,2:28)]
setwd(res.dir)
write_xlsx(M3.screen.df, path="M3_PairwiseMRD_Screen_2023-07-16.xlsx")

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
M3.screen.df2 <- M3.screen.df2[,c(1,29,2:28)]
setwd(res.dir)
write_xlsx(M3.screen.df2, path="M3_PairwiseMRD_Screen_Factor_2023-07-16.xlsx")

## M3 genetic subtype
library(janitor)
M3.genetic.subtype2 <- clean_names(M3.genetic.subtype)
rownames(M3.genetic.subtype2) <- rownames(M3.genetic.subtype)

M3g.screen.df <- uni_screen(annot, M3.genetic.subtype2)
M3g.screen.df$MolecularFeatureOld <- colnames(M3.genetic.subtype)
M3g.screen.df <- M3g.screen.df[,c(1, 29, 2:28)]
setwd(res.dir)
write_xlsx(M3g.screen.df, path="M3_genetic_PairwiseMRD_Screen_2023-07-16.xlsx")

## M3 subsubtype
library(janitor)
M3.subsubtype2 <- clean_names(M3.subsubtype)
rownames(M3.subsubtype2) <- rownames(M3.subsubtype)

M3s.screen.df <- uni_screen(annot, M3.subsubtype2)
M3s.screen.df$MolecularFeatureOld <- colnames(M3.subsubtype)
M3s.screen.df <- M3s.screen.df[,c(1,29,2:28)]
setwd(res.dir)
write_xlsx(M3s.screen.df, path="M3_subsubtype_PairwiseMRD_Screen_2023-07-16.xlsx")

### M4
# Remove problematic characters from column names
M4.pathway2 <- clean_names(M4.pathway)
M4.data <- M4.pathway2
rownames(M4.data) <- rownames(M4.pathway)

M4.screen.df <- uni_screen(annot, M4.data)
M4.screen.df$MolecularFeatureOld <- colnames(M4.pathway)
M4.screen.df <- M4.screen.df[,c(1,29,2:28)]
setwd(res.dir)
write_xlsx(M4.screen.df, path="M4_PairwiseMRD_Screen_2023-07-16.xlsx")

### M7 IP
M7.IP2 <- clean_names(M7.IP)
rownames(M7.IP2) <- rownames(M7.IP)

M7.screen.df <- uni_screen(annot, M7.IP2)
M7.screen.df$MolecularFeatureOld <- colnames(M7.IP)
M7.screen.df <- M7.screen.df[,c(1,29,2:28)]
setwd(res.dir)
library(writexl)
write_xlsx(M7.screen.df, path="M7_PairwiseMRD_Screen_2023-07-16.xlsx")


### M5 Genes
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


M5.screen.df <- uni_screen_gt(annot, M5.AllLesions.variants2, gene.list)
res <- M5.screen.df$Results
gene.names2 <- gene.names
colnames(gene.names2)[2] <- "Gene"
res2 <- dplyr::full_join(res, gene.names2, by="Gene")
res3 <- res2[,c(1,23,2:22)]
M5.screen.df <- res3
setwd(res.dir)
write.csv(M5.screen.df, file="M5_genes_PairwiseMRD_Screen_2023-07-16.csv")

### M5 Variants
mrd.NegRD.sig.genes <- res3[which(res3$MRD.NegRD.P<0.1),]$Gene
mrd.NegIF.sig.genes <- res3[which(res3$MRD.NegIF.P<0.1),]$Gene
mrd.RDIF.sig.genes <- res3[which(res3$MRD.RDIF.P<0.1),]$Gene
sig.genes.comb <- c(mrd.NegRD.sig.genes, mrd.NegIF.sig.genes, mrd.RDIF.sig.genes)
sig.genes.un <- unique(sig.genes.comb)
sig.gene.list <- gene.list[names(gene.list) %in% sig.genes.un]
sig.variants <- unlist(sig.gene.list)
# Run univariate screen on filtered variant list
M5.vars.sig <- M5.AllLesions.variants2[,which(colnames(M5.AllLesions.variants2) %in% sig.variants)]
M5.vars.sig.oldnames <- colnames(M5.AllLesions.variants[,which(colnames(M5.AllLesions.variants2) %in% sig.variants)])

M5.screen.df2 <- uni_screen(annot, M5.vars.sig)
M5.screen.df2$MolecularFeaturesOld <- M5.vars.sig.oldnames
M5.screen.df2 <- M5.screen.df2[,c(1,29,2:28)]
setwd(res.dir)
write.csv(M5.screen.df2, file="M5_variants_PairwiseMRD_Screen_2023-07-16.csv")

### Gene Expression - Run on HPC

# iqr.scale.func <- function(gene.dat){
#   m <- median(as.numeric(gene.dat))
#   s <- IQR(gene.dat)
#   scale.gene.dat <- (gene.dat-m)/s
#   return(scale.gene.dat)
# }
# gexp.scale <- apply(gexp, 1, iqr.scale.func)
# nan.vec <- colSums(is.nan(gexp.scale))
# length(which(nan.vec>0)) #1019
# gexp.scale.filt <- gexp.scale[,-which(nan.vec>0)]
# all.equal(rownames(gexp.scale.filt), rownames(annot))
# 
# gexp.screen.df <- uni_screen(annot, gexp.scale.filt)
# setwd(res.dir)
# write.csv(gexp.screen.df, file="Gexp_Uni_Screen_2023-06-22.csv")
