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
#Path_file="/public/workspace/liuqzh/gastric_cancer/GSE150290_inferCNV_NMF_result/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE150290_sample_data/"

#Anuja
#Path_file="/public/workspace/liuqzh/gastric_cancer/Anuja_inferCNV_NMF_result/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/Anuja_sample_data/"

#CRA002586
#Path_file="/public/workspace/liuqzh/gastric_cancer/CRA002586_inferCNV_NMF_result/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/CRA002586_sample_data/"

#GSE167297
#Path_file="/public/workspace/liuqzh/gastric_cancer/GSE167297_inferCNV_NMF_result/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE167297_sample_data/"

#GSE183904
#Path_file="/public/workspace/liuqzh/gastric_cancer/GSE183904_inferCNV_NMF_result/"
#Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE183904_sample_data/"

#####GSE150290
{
Path_file="/public/workspace/liuqzh/gastric_cancer/GSE150290_inferCNV_NMF_result/"
Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE150290_sample_data/"

cell_name=NULL
for(Sample in list.files(Path_file)){
	print(Sample)
	setwd(paste(Path_file,Sample,sep=''))
	dat<-readRDS(paste(Celltype_file,Sample,'/Single_cell_celltype_CNA_new.rds',sep=''))	
	dat<-subset(dat,subset=CNA_celltype %in% 'Tumor_cell')
	
	tmp_data=cbind(dat@meta.data$cell_idents,dat@meta.data$CNA_celltype)
	cell_name=rbind(cell_name,tmp_data)
}
colnames(cell_name)=c('cell_idents','CNA_celltype')
GSE150290=cell_name
}

#####Anuja
####
{
Path_file="/public/workspace/liuqzh/gastric_cancer/Anuja_inferCNV_NMF_result/"
Celltype_file="/public/workspace/liuqzh/gastric_cancer/Anuja_sample_data/"

cell_name=NULL
for(Sample in list.files(Path_file)){
	print(Sample)
	setwd(paste(Path_file,Sample,sep=''))
	dat<-readRDS(paste(Celltype_file,Sample,'/Single_cell_celltype_CNA_new.rds',sep=''))	
	dat<-subset(dat,subset=CNA_celltype %in% 'Tumor_cell')
	
	tmp_data=cbind(dat@meta.data$cell_idents,dat@meta.data$CNA_celltype)
	cell_name=rbind(cell_name,tmp_data)
}
colnames(cell_name)=c('cell_idents','CNA_celltype')
Anuja=cell_name
}

#####CRA002586
####
{
Path_file="/public/workspace/liuqzh/gastric_cancer/CRA002586_inferCNV_NMF_result/"
Celltype_file="/public/workspace/liuqzh/gastric_cancer/CRA002586_sample_data/"

cell_name=NULL
for(Sample in list.files(Path_file)){
	print(Sample)
	setwd(paste(Path_file,Sample,sep=''))
	dat<-readRDS(paste(Celltype_file,Sample,'/Single_cell_celltype_CNA_new.rds',sep=''))	
	dat<-subset(dat,subset=CNA_celltype %in% 'Tumor_cell')
	
	tmp_data=cbind(dat@meta.data$cell_idents,dat@meta.data$CNA_celltype)
	cell_name=rbind(cell_name,tmp_data)
}
colnames(cell_name)=c('cell_idents','CNA_celltype')
CRA002586=cell_name
}

#####GSE167297
####
{
Path_file="/public/workspace/liuqzh/gastric_cancer/GSE167297_inferCNV_NMF_result/"
Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE167297_sample_data/"

cell_name=NULL
for(Sample in list.files(Path_file)){
	print(Sample)
	setwd(paste(Path_file,Sample,sep=''))
	dat<-readRDS(paste(Celltype_file,Sample,'/Single_cell_celltype_CNA_new.rds',sep=''))	
	dat<-subset(dat,subset=CNA_celltype %in% 'Tumor_cell')
	
	tmp_data=cbind(dat@meta.data$cell_idents,dat@meta.data$CNA_celltype)
	cell_name=rbind(cell_name,tmp_data)
}
colnames(cell_name)=c('cell_idents','CNA_celltype')
GSE167297=cell_name
}

#####GSE183904
####
{
Path_file="/public/workspace/liuqzh/gastric_cancer/GSE183904_inferCNV_NMF_result/"
Celltype_file="/public/workspace/liuqzh/gastric_cancer/GSE183904_sample_data/"

cell_name=NULL
for(Sample in list.files(Path_file)){
	print(Sample)
	setwd(paste(Path_file,Sample,sep=''))
	dat<-readRDS(paste(Celltype_file,Sample,'/Single_cell_celltype_CNA_new.rds',sep=''))	
	dat<-subset(dat,subset=CNA_celltype %in% 'Tumor_cell')
	
	tmp_data=cbind(dat@meta.data$cell_idents,dat@meta.data$CNA_celltype)
	cell_name=rbind(cell_name,tmp_data)
}
colnames(cell_name)=c('cell_idents','CNA_celltype')
GSE183904=cell_name
}


data_cell_names=rbind(GSE150290,Anuja,CRA002586,GSE167297,GSE183904)
#saveRDS(data_cell_names,'/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_data_cell_names_CNA.rds')






