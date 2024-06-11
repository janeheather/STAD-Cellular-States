bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE210249/')
GSE210249<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE210249/GSE210249_data.txt',header=T,row.names=1,sep='\t')

eff_length <- read.table("/public/workspace/yumiao/PDAC/lnrna/efflen_symbol.txt",header=T,sep = "\t")

GSE210249$geneid_efflen.geneid=rownames(GSE210249)
##data为count矩阵
data <- merge(GSE210249,eff_length,by='geneid_efflen.geneid')

rownames(data)=data$exon_symbol.gene_name
GSE210249_data=data[,-c(1,8,9)]
GSE210249_data <- log2(GSE210249_data+1)

saveRDS(GSE210249_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE210249/GSE210249_data.rds')

