
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(dplyr)
library(ggplot2)
library(ggtern)
library(GSVA)
library(org.Mm.eg.db)
library(clusterProfiler)

set.seed(520)
GSE184613_data=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE184613/GSE184613_data.txt",sep = "\t",header = T)

eg <- bitr(GSE184613_data[,1],fromType = 'ENSEMBL',
           toType = c('SYMBOL','ENTREZID'),
           OrgDb='org.Mm.eg.db',
)
colnames(GSE184613_data)[1]='ENSEMBL'
GSE184613_data=as.data.frame(GSE184613_data)
eg=as.data.frame(eg)

GSE184613_data<-merge(eg,GSE184613_data,by='ENSEMBL')
##计算count转tpm
x <- GSE184613_data[,5:ncol(GSE184613_data)] /(GSE184613_data$Length/1000)

tpm = t( t(x) / colSums(x) ) * 1e6 
tpm_log <- as.data.frame(log2(tpm+1))
tpm_log$SYMBOL=GSE184613_data$SYMBOL

final<-aggregate(tpm_log[,1:6],list(tpm_log[,7]),mean)
rownames(final)=final[,1]
final=final[,-1]

rownames(final)=toupper(rownames(final))
GSE184613_data=final
#saveRDS(GSE184613_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE184613/GSE184613_data.rds')

