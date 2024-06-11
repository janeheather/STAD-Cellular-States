bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ComplexHeatmap)

dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(dat) <- "RNA"
dat=subset(dat,cells=which(dat$celltype=='Fibroblast'|dat$celltype=='B cell'|
						   dat$celltype=='Endothelial'|dat$celltype=='Epithelium'|
						   dat$celltype=='Macrophage'|dat$celltype=='NK cell'|
						   dat$celltype=='T cell'|dat$celltype=='Mast cell'))

scdata_score<-readRDS('/public/workspace/liuqzh/gastric_cancer/scalop-master/scdata_score.rds')
scdata_score=as.matrix(scdata_score)

dat_result=scdata_score[intersect(rownames(dat@meta.data),rownames(scdata_score)),]
dat@meta.data=cbind(dat@meta.data,dat_result)

FeaturePlot(dat,features=c('Mac_MP_1'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Mac_MP_2'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Mac_MP_3'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))


FeaturePlot(dat,features=c('T_cell_MP_1'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('T_cell_MP_2'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('T_cell_MP_3'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('T_cell_MP_4'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('T_cell_MP_5'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('T_cell_MP_6'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))

FeaturePlot(dat,features=c('Mast_MP_1'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Mast_MP_2'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Mast_MP_3'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Mast_MP_4'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))

FeaturePlot(dat,features=c('B_MP_1'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('B_MP_2'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('B_MP_3'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('B_MP_4'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('B_MP_5'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('B_MP_6'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))

FeaturePlot(dat,features=c('Endo_MP_1'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Endo_MP_2'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Endo_MP_3'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Endo_MP_4'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Endo_MP_5'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Endo_MP_6'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))


DimPlot(dat,group.by="subtype",reduction = "umap",cols=colors,label=TRUE,raster=FALSE)
FeaturePlot(dat,features=c('Mac_MP_1'),reduction='tsne',coord.fixed=TRUE,order=TRUE,pt.size=0.1,cols = c("#a1a3a6", "#f47920"))
###############

bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ComplexHeatmap)

dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(dat) <- "RNA"
dat=subset(dat,cells=which(dat$celltype=='Fibroblast'|dat$celltype=='B cell'|
						   dat$celltype=='Endothelial'|dat$celltype=='Epithelium'|
						   dat$celltype=='Macrophage'|dat$celltype=='NK cell'|
						   dat$celltype=='T cell'|dat$celltype=='Mast cell'))

scdata_score<-readRDS('/public/workspace/liuqzh/gastric_cancer/scalop-master/scdata_score.rds')
scdata_score=as.matrix(scdata_score)

dat_result=scdata_score[intersect(rownames(dat@meta.data),rownames(scdata_score)),]
dat@meta.data=cbind(dat@meta.data,dat_result)

table(dat@meta.data$celltype)
B cell     Endothelial  Epithelium  Fibroblast  Macrophage   Mast cell
  52266        9818       32120       17541       21271        5313
NK cell      T cell
  22125       50117

######B cell
dat_celltype_tmp=dat@meta.data[which(dat@meta.data$celltype == 'B cell'),]
dat_celltype_tmp=dat_celltype_tmp[,grep('^B_MP',colnames(dat_celltype_tmp))]
NMF_class_B_cell=as.data.frame(paste('B_MP_',apply(dat_celltype_tmp,1,which.max),sep=''))
rownames(NMF_class_B_cell)=rownames(dat_celltype_tmp)
colnames(NMF_class_B_cell)='NMF_class'

######Endothelial
dat_celltype_tmp=dat@meta.data[which(dat@meta.data$celltype == 'Endothelial'),]
dat_celltype_tmp=dat_celltype_tmp[,grep('^Endo',colnames(dat_celltype_tmp))]
NMF_class_Endothelial=as.data.frame(paste('Endo_MP_',apply(dat_celltype_tmp,1,which.max),sep=''))
rownames(NMF_class_Endothelial)=rownames(dat_celltype_tmp)
colnames(NMF_class_Endothelial)='NMF_class'

######Epithelium
dat_celltype_tmp=dat@meta.data[which(dat@meta.data$celltype == 'Epithelium'),]
dat_celltype_tmp=dat_celltype_tmp[,grep('^Maligant_MP',colnames(dat_celltype_tmp))]
NMF_class_Mal_cell=as.data.frame(paste('Maligant_MP_',apply(dat_celltype_tmp,1,which.max),sep=''))
rownames(NMF_class_Mal_cell)=rownames(dat_celltype_tmp)
colnames(NMF_class_Mal_cell)='NMF_class'

######Fibroblast
dat_celltype_tmp=dat@meta.data[which(dat@meta.data$celltype == 'Fibroblast'),]
dat_celltype_tmp=dat_celltype_tmp[,grep('^Fib_MP',colnames(dat_celltype_tmp))]
NMF_class_Fib_cell=as.data.frame(paste('Fib_MP_',apply(dat_celltype_tmp,1,which.max),sep=''))
rownames(NMF_class_Fib_cell)=rownames(dat_celltype_tmp)
colnames(NMF_class_Fib_cell)='NMF_class'

######Macrophage
dat_celltype_tmp=dat@meta.data[which(dat@meta.data$celltype == 'Macrophage'),]
dat_celltype_tmp=dat_celltype_tmp[,grep('^Mac_MP',colnames(dat_celltype_tmp))]
NMF_class_Mac_cell=as.data.frame(paste('Mac_MP_',apply(dat_celltype_tmp,1,which.max),sep=''))
rownames(NMF_class_Mac_cell)=rownames(dat_celltype_tmp)
colnames(NMF_class_Mac_cell)='NMF_class'

######Mast cell
dat_celltype_tmp=dat@meta.data[which(dat@meta.data$celltype == 'Mast cell'),]
dat_celltype_tmp=dat_celltype_tmp[,grep('^Mast_MP',colnames(dat_celltype_tmp))]
NMF_class_Mast_cell=as.data.frame(paste('Mast_MP_',apply(dat_celltype_tmp,1,which.max),sep=''))
rownames(NMF_class_Mast_cell)=rownames(dat_celltype_tmp)
colnames(NMF_class_Mast_cell)='NMF_class'

######NK cell/T cell
dat_celltype_tmp=dat@meta.data[which(dat@meta.data$celltype %in% c('NK cell','T cell')),]
dat_celltype_tmp=dat_celltype_tmp[,grep('^T_cell_MP',colnames(dat_celltype_tmp))]
NMF_class_T_cell=as.data.frame(paste('T_cell_MP_',apply(dat_celltype_tmp,1,which.max),sep=''))
rownames(NMF_class_T_cell)=rownames(dat_celltype_tmp)
colnames(NMF_class_T_cell)='NMF_class'

NMF_class=rbind(NMF_class_B_cell,NMF_class_Endothelial,NMF_class_Mal_cell,
				NMF_class_Fib_cell,NMF_class_Mac_cell,NMF_class_Mast_cell,
				NMF_class_T_cell)

NMF_class=NMF_class[intersect(rownames(dat@meta.data),rownames(NMF_class)),]
dat@meta.data=cbind(dat@meta.data,NMF_class)

dat@meta.data$NMF_subtype=''
dat@meta.data$NMF_subtype[dat@meta.data$NMF_class %in% c('Mac_MP_3','Endo_MP_1','Endo_MP_7',
														 'Endo_MP_5','Fib_MP_1','Fib_MP_4',
														 'Maligant_MP_2','Mast_MP_4','B_MP_6')]='Com_1'
														 
dat@meta.data$NMF_subtype[dat@meta.data$NMF_class %in% c('B_MP_5','T_cell_MP_1','Mast_MP_3',
														 'Endo_MP_6','Maligant_MP_1','T_cell_MP_3',
														 'B_MP_4','T_cell_MP_4','B_MP_3',
														 'Mac_MP_2','Endo_MP_3')]='Com_2'
														 
dat@meta.data$NMF_subtype[dat@meta.data$NMF_class %in% c('Mac_MP_1','Fib_MP_2')]='Com_3'

dat@meta.data$NMF_subtype[dat@meta.data$NMF_class %in% c('Endo_MP_2','Mast_MP_2','Fib_MP_3',
														 'T_cell_MP_5','Maligant_MP_3','Mast_MP_1',
														 'T_cell_MP_6','B_MP_2','T_cell_MP_2',
														 'T_cell_MP_2','B_MP_1','Endo_MP_4')]='Com_4'


saveRDS(dat,'/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_NMF_class.rds')

FeaturePlot(dat,features=c('Fib_MP_1'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Fib_MP_2'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Fib_MP_3'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))
FeaturePlot(dat,features=c('Fib_MP_4'),reduction='umap',coord.fixed=TRUE,min.cutoff=2,max.cutoff=5,order=TRUE,pt.size=0.5,cols = c("#a1a3a6", "#f47920"))























