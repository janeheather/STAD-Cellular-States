bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE246963/')
GSE246963<-read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE246963/GSE246963_data.csv',header=T,sep=',')
eff_length <- read.table("/public/workspace/yumiao/PDAC/lnrna/efflen_symbol.txt",header=T,sep = "\t")
colnames(GSE246963)[1]='geneid_efflen.geneid'

##data为count矩阵
GSE246963 <- merge(eff_length,GSE246963,by='geneid_efflen.geneid')
final<-aggregate(GSE246963[,3:ncol(GSE246963)],list(GSE246963[,2]),mean)
rownames(final)=final[,1]

GSE246963_data=final[,-c(1:2)]

saveRDS(GSE246963_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE246963/GSE246963_data.rds')

