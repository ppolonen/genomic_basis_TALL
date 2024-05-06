rm(list=ls());              # start with clean slate in R
options(stringsAsFactors=F) # turn off the default in R

knitr::opts_chunk$set(echo=T,error=T,eval=T)

# Specify a directory on your local machine
setwd("~/GRIN_results")

library(readxl)
library(openxlsx)
library(writexl)
library(circlize)
library(biomaRt)
library(ComplexHeatmap)
library(ggplot2)
library(forcats)
library(EnvStats)
library(tidyverse)
library(lubridate)
library(GSDA)
library(stringr)
library(survival)
library(dplyr)

# GRIN2 source code:
source("https://raw.githubusercontent.com/stjude/TALL-example/main/GRIN2.0.ALEX.library.09.29.2022.R")

## 1) Obtain Genomic Annotations
hg38.ann=get.ensembl.annotation("Human_GRCh38") 
# "Human_GRCh38" can be used to retrieve data for hg38
hg38.gene.annotation=hg38.ann$gene.annotation
# Gene annotation data that include around 20,000 coding genes and 25,000 Non-coding processed 
# transcripts such as lncRNAs, miRNAs, snRNA and snoRNAs
hg38.reg.annotation=hg38.ann$reg.annotation
# Annotation data for regulatory features retrieved from ensembl regulatory build that include 
# around 600,000 feauters (promoters, enhancer, TF and CTCF binding sites, etc...)
# Ensembl imports publicly available data from different large epigenomic consortia that includes 
# ENCODE, Roadmap Epigenomics and Blueprint (118 epigenome)
symbol_ensembl=cbind.data.frame(gene.name=hg38.gene.annotation$gene.name,
                                gene=hg38.gene.annotation$gene)

hg38.chrom.size=get.chrom.length("Human_GRCh38")


##*****************************************************************************************
## Running notes: Data can be downloaded from synapse syn54032669
##*****************************************************************************************

# load coding lesion data:
load("All_GRIN_lesions_variantID_coding.Rdata")

GRIN.results=grin.stats(lesions[,1:6],
                        hg38.gene.annotation, 
                        hg38.chrom.size)

save(list=c("GRIN.results"), file = "GRIN_results_exonic.Rdata")

# load regulatory lesion data, see methods for data processing:
load("All_GRIN_lesions_variantID_regulatory.Rdata")
GRIN.results.regulatory=grin.stats(lesions[,1:5],
                                   hg38.reg.annotation, 
                                   hg38.chrom.size)

save(list=c("GRIN.results.regulatory"), file = "GRIN_results_regulatory_allLesions.Rdata")
