bytlib load R-4.0.2
bytlib load gcc
R

library(NMF)
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(scran)
library(patchwork)
library(ggplot2)

######
#GSE150290
#Path_file="/public/workspace/liuqzh/gastric_cancer/GSE150290_celltype_endo_NMF/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE150290_sample_data/"

#Anuja
#Path_file="/public/workspace/liuqzh/gastric_cancer/Anuja_celltype_endo_NMF/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/Anuja_sample_data/"

#CRA002586
#Path_file="/public/workspace/liuqzh/gastric_cancer/CRA002586_celltype_endo_NMF/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/CRA002586_sample_data/"

#GSE167297
#Path_file="/public/workspace/liuqzh/gastric_cancer/GSE167297_celltype_endo_NMF/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE167297_sample_data/"

#GSE183904
#Path_file="/public/workspace/liuqzh/gastric_cancer/GSE183904_celltype_endo_NMF/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE183904_sample_data/"

Path_file="/public/workspace/liuqzh/gastric_cancer/GSE183904_celltype_endo_NMF/"
Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE183904_sample_data/"

for(Sample in list.files(Path_file)){
	print(Sample)
	setwd(paste(Path_file,Sample,sep=''))
	dat<-readRDS(paste(Celltype_file,Sample,'/Single_cell_celltype.rds',sep=''))
	dat<-subset(dat,subset=single_cell_celltype %in% 'Endothelial')
	##########编码蛋白处理
	gene.type=read.table("/public/workspace/liuqzh/lqz_ref/refdata-cellranger-hg19-3.0.0/genes/gene.type",header=F)
	colnames(gene.type)=c("ENSG_id","Gene_id","Gene_type")

	gene.type=gene.type[which(gene.type$Gene_type=="protein_coding"),]
	gene.type=gene.type[which(!duplicated(gene.type$Gene_id)),]
	###
	dat_sample<-dat
	data=as.matrix(dat_sample@assays$RNA@data)

	rownames(gene.type)=gene.type$Gene_id
	data=data[intersect(rownames(data),rownames(gene.type)),]

	ratio=apply(data,1,function(x){sum(x>0)/length(x)})
	dest_ratio=as.data.frame(ratio)
	p_ratio=ggplot(dest_ratio,aes(x=ratio))+geom_density()

	express=apply(data,1,function(x){mean(x[x>0])})
	dest_exp=as.data.frame(express)
	p_exp=ggplot(dest_exp,aes(x=express))+geom_density()

	p_ratio+p_exp

	dest_exp=dest_exp[which(is.nan(dest_exp[,1])==F),,drop=F]
	dest_exp=arrange(dest_exp,desc(express))
	express=dest_exp[c(1:ceiling(length(dest_exp[which(is.nan(dest_exp[,1])==F),])*0.4)),,drop=F]###前40%

	####参数选择####
	data=data[which(rownames(data) %in% names(which(0.75 > ratio & ratio > 0.05))),]
	data=data[which(rownames(data) %in% rownames(express)),]
	################

	data_RNA=data
	data_count=as.matrix(dat_sample@assays$RNA@counts)

	data=data_count[intersect(rownames(data_RNA),rownames(data_count)),]

	#MAD=apply(data,1,function(x){mad(x[x>0])})
	write.table(t(data),paste(Path_file,Sample,'/','0.4_0.75_0.05_count_data.txt',sep=''),col.names=T, sep='\t')
}
