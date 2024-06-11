#############################
bytlib load R-4.0.2
bytlib load gcc
R

library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)
library(Seurat)
library(dplyr)

STAD_P<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/STAD_P_E_T_N.rds')
DefaultAssay(STAD_P)<-'RNA'
dat=STAD_P

TPM=as.data.frame(dat[["RNA"]]@data)
TPM<-TPM[rowSums(TPM)>0,]
TPM<-as.matrix(TPM)

genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/subtype_III.csv",sep = ",",row.names = 1,header = F)
geneset<-c()
for (i in 1:nrow(genelist)) {
	aa<-genelist[i,]
	aa<-aa[!is.na(aa)]
	aa<-aa[!duplicated(aa)]
	aa<-aa[aa%in%rownames(TPM)]
	aa<-list(as.character(as.matrix(aa)))
	geneset<-c(geneset,aa)
}
names(geneset)<-rownames(genelist)
colCenter = function(m, by = 'mean') {
  m = as.matrix(m)
  if (by == 'mean')  by = T
  else if (by == 'median') by = matrixStats::colMedians(m)
  else stopifnot(is.numeric(by) & length(by) == ncol(m))
  scale(m, center = by, scale = F)
}
sc = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)
sc=as.data.frame(sc)
#sc$PLF=sc$PLF1+sc$PLF2
result=NULL
data_meta=as.matrix(sc)
#data_meta=data_meta[,3:5]

Subtype=NULL
for(i in 1:nrow(data_meta)){
	ch_max=pmax(data_meta[i,1],data_meta[i,2],data_meta[i,3])
	if(ch_max == data_meta[i,1]){
		Subtype=c(Subtype,'MP_1')
	}	
	if(ch_max == data_meta[i,2]){
		Subtype=c(Subtype,'MP_2')
	}	
	if(ch_max == data_meta[i,3]){
		Subtype=c(Subtype,'MP_3')
	}	
}
data_meta=as.data.frame(data_meta)	
data_meta$Subtype=Subtype

dat@meta.data=cbind(dat@meta.data,data_meta)

DimPlot(dat,group.by="Subtype",reduction='umap',pt.size=0.1,label = TRUE)

pdf('/public/workspace/liuqzh/gastric_cancer/figure2/MP_plot.pdf')
FeaturePlot(dat,features=c('MP_1'),reduction='umap',coord.fixed=TRUE,order=TRUE,pt.size=0.1,min.cutoff=0,cols = c("#a1a3a6", "#BD3934"))
FeaturePlot(dat,features=c('MP_2'),reduction='umap',coord.fixed=TRUE,order=TRUE,pt.size=0.1,min.cutoff=0,cols = c("#a1a3a6", "#BD3934"))
FeaturePlot(dat,features=c('MP_3'),reduction='umap',coord.fixed=TRUE,order=TRUE,pt.size=0.1,min.cutoff=0,cols = c("#a1a3a6", "#BD3934"))
dev.off()

saveRDS(dat, '/public/workspace/liuqzh/gastric_cancer/figure2/STAD_Epithelium_Subtype_Score.rds')

pdf('/public/workspace/liuqzh/gastric_cancer/single_cell/EMT.pdf',height=6,width=6)
FeaturePlot(dat,features=c('EMT'),reduction='umap',min.cutoff=1,max.cutoff=2,order=TRUE,pt.size=0.3,cols = c("#a1a3a6", "#f47920"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/single_cell/PLF.pdf',height=6,width=6)
FeaturePlot(dat,features=c('PLF'),reduction='umap',min.cutoff=1,max.cutoff=2,order=TRUE,pt.size=0.3,cols = c("#a1a3a6", "#f47920"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/single_cell/NLT.pdf',height=6,width=6)
FeaturePlot(dat,features=c('NLT'),reduction='umap',min.cutoff=0.7,max.cutoff=1.4,order=TRUE,pt.size=0.3,cols = c("#a1a3a6", "#f47920"))
dev.off()


