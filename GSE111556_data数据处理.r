
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(dplyr)
library(ggplot2)
library(ggtern)
library(GSVA)

set.seed(520)
GSE111556_data=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE111556/GSE111556_data.txt",sep = "\t",header = T)
final<-aggregate(GSE111556_data[,2:ncol(GSE111556_data)],list(GSE111556_data[,1]),mean)
GSE111556_data=final
rownames(GSE111556_data)=GSE111556_data[,1]

eff_length <- read.table("/public/workspace/yumiao/PDAC/lnrna/efflen_symbol.txt",header=T,sep = "\t")
colnames(GSE111556_data)[1]='exon_symbol.gene_name'
eff_length<- eff_length[!duplicated(eff_length$exon_symbol.gene_name),]

GSE111556_data<-merge(eff_length,GSE111556_data,by='exon_symbol.gene_name')
rownames(GSE111556_data)<-GSE111556_data$exon_symbol.gene_name

##计算count转tpm
x <- GSE111556_data[,4:ncol(GSE111556_data)] /(GSE111556_data$geneid_efflen.efflen/1000)

tpm = t( t(x) / colSums(x) ) * 1e6 

tpm_log <- log2(tpm+1)
GSE111556_data=tpm_log

#saveRDS(GSE111556_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE111556/GSE111556_data.rds')
GSE111556_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE111556/GSE111556_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE111556_data<-GSE111556_data[rowSums(GSE111556_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE111556_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE111556_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE111556_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')

