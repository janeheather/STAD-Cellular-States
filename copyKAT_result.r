
bytlib load R-4.0.2
bytlib load gcc
R

library(NMF)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(Seurat)
library(copykat)
library(dplyr)

STAD_Epithelium_Subtype_Score<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')

counts = as.matrix(STAD_Epithelium_Subtype_Score@assays$RNA@counts)
copyKAT_results = copykat(rawmat = counts,ngene.chr = 5,sam.name = 'STAD',n.cores = 40)
# 整合预测结果到Seurat对象中
STAD_Epithelium_Subtype_Score$cell.names<-rownames(STAD_Epithelium_Subtype_Score@meta.data)
STAD_Epithelium_Subtype_Score@meta.data<-merge(copyKAT_results$prediction,STAD_Epithelium_Subtype_Score@meta.data,by='cell.names')
rownames(STAD_Epithelium_Subtype_Score@meta.data)<-STAD_Epithelium_Subtype_Score@meta.data$cell.names
p1=DimPlot(STAD_Epithelium_Subtype_Score,group.by="celltype",reduction = "tsne",label=TRUE,raster=FALSE)
p2=DimPlot(STAD_Epithelium_Subtype_Score,group.by="copykat.pred",reduction = "tsne",label=TRUE,raster=FALSE)
p1+p2
saveRDS(STAD_Epithelium_Subtype_Score,"/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype_copyKAT_results.rds")

data_cell_names=readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_data_cell_names_CNA.rds')
rownames(data_cell_names)=data_cell_names[,1]
dat<-subset(STAD_Epithelium_Subtype_Score,subset=cell.names %in% rownames(data_cell_names))
p3=DimPlot(dat,group.by="celltype",reduction = "tsne",label=TRUE,raster=FALSE)
p4=DimPlot(dat,group.by="copykat.pred",reduction = "tsne",label=TRUE,raster=FALSE)
p3+p4



