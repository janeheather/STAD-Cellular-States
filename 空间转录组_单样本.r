module load hdf5-1.8.13
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(hdf5r)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)

setwd("/public/workspace/liuqzh/gastric_cancer/sp_result/")
expr <- "/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE251950_SP/P21_00732_LI_SING/21_00732_LI_SING_filtered_feature_bc_matrix.h5"
expr.mydata <- Seurat::Read10X_h5(filename =  expr )
mydata <- Seurat::CreateSeuratObject(counts = expr.mydata, project = 'STAD', assay = 'Spatial')
mydata$slice <- 1
mydata$region <- 'P21_00732_LI_SING' #命名

plot1 <- VlnPlot(mydata, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(mydata, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)

mydata <- SCTransform(mydata, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
# also run standard log normalization for comparison
mydata<- NormalizeData(mydata, verbose = FALSE, assay = "Spatial")

# Computes the correlation of the log normalized data and sctransform residuals with the number
# of UMIs
mydata <- GroupCorrelation(mydata, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
mydata <- GroupCorrelation(mydata, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
p1 <- GroupCorrelationPlot(mydata, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") + 
    theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(mydata, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") + 
    theme(plot.title = element_text(hjust = 0.5))
plot_grid(p1, p2)

mydata <- RunPCA(mydata, assay = "SCT", verbose = FALSE)
mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:30,)
mydata <- FindClusters(mydata, verbose = FALSE,resolution = 1)
mydata <- RunTSNE(mydata, reduction = "pca", dims = 1:30)


pdf('/public/workspace/yumiao/scrabble_score.pdf',width=6,height=6,useDingbats=F)
DimPlot(mydata, reduction = "tsne", label = TRUE)
dev.off()

FeaturePlot(mydata,features=c('PTPRC','EPCAM','COL1A2','VWF','CHGA'),reduction='tsne',coord.fixed=TRUE,order=TRUE,pt.size=0.3,cols = c("#a1a3a6", "#f47920"))



FeaturePlot(mydata,features=c('EPCAM'),reduction='tsne',coord.fixed=TRUE,order=TRUE,pt.size=0.3,cols = c("#a1a3a6", "#f47920"))

p2 <- SpatialDimPlot(mydata, label = TRUE, label.size = 3)
plot_grid(p1, p2)
dev.off()