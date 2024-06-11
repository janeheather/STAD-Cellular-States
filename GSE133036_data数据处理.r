bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE133036/')
GSE133036<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE133036/GSE133036_data.txt',header=T,sep='\t')

final<-aggregate(GSE133036[,2:ncol(GSE133036)],list(GSE133036[,1]),mean)
rownames(final)<-final[,1]
GSE133036=final[,-1]

eff_length <- read.table("/public/workspace/yumiao/PDAC/lnrna/efflen_symbol.txt",header=T,sep = "\t")

GSE133036$exon_symbol.gene_name=rownames(GSE133036)
##data为count矩阵
data <- merge(eff_length,GSE133036,by='exon_symbol.gene_name')

##基因长度文件symbol id去重
data<-aggregate(data[,4:ncol(data)],list(data[,2]),max)
eff_length<- eff_length[!duplicated(eff_length$exon_symbol.gene_name),]
colnames(data)[1]<-'geneid_efflen.geneid'
data<-merge(eff_length,data,by='geneid_efflen.geneid')
rownames(data)<-data$exon_symbol.gene_name

##计算count转tpm
x <- data[,4:ncol(data)] /(data$geneid_efflen.efflen/1000)

tpm = t( t(x) / colSums(x) ) * 1e6 

tpm_log <- log2(tpm+1)
GSE133036_data=tpm_log
saveRDS(GSE133036_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE133036/GSE133036_data.rds')

