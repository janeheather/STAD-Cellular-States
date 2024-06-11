############CytoTRACE#########
{
bytlib load R-4.0.2
bytlib load gcc
R

library(Seurat)
library(CytoTRACE)
library(monocle)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)
set.seed(520)
#读取scRNA
dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(dat) <- "RNA"

setwd('/public/workspace/yumiao/TE_GAST/')
data_single_cell=dat
DimPlot(data_single_cell, reduction = "tsne",label=TRUE,raster=FALSE)
DimPlot(data_single_cell, reduction = "umap",label=TRUE,raster=FALSE)

exprdat <- as.matrix(GetAssayData(data_single_cell[["RNA"]], slot='counts'))
results_CytoTRACE <- CytoTRACE(exprdat)

phenotype=as.matrix(data_single_cell@meta.data)
phenotype=t(phenotype)
phenotype=phenotype[16,]

emb=as.data.frame(Embeddings(data_single_cell,reduction = "umap"))
p=plotCytoTRACE(results_CytoTRACE,emb = emb,phenotype=phenotype)

data_single_cell_metadata=data_single_cell@meta.data
data_single_cell_metadata$CytoTRACE_Score<-as.numeric(results_CytoTRACE$CytoTRACE)

test<-data_single_cell_metadata    
test$cluster<-factor(test$Subtype)

test=test[which(test$Subtype=="MP_1"|test$Subtype=="MP_2"|test$Subtype=="MP_3"),]
test$cluster<-factor(test$Subtype,levels=c("MP_1","MP_2","MP_3"))
Class <- list(c("MP_1","MP_2"),c("MP_2","MP_3"),c("MP_1","MP_3"))

colors=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")

p1 =ggplot(test,aes(x=cluster,y=CytoTRACE_Score)) +
	geom_violin(aes(fill=Subtype),alpha=0.7,outlier.alpha=0,width=0.7)+	
	geom_boxplot(width=0.2,alpha=0.7,outlier.alpha=0)+   
	labs(x = '')+
	scale_y_continuous(name = "Score")+
	ggtitle("CytoTRACE_Score")+
	scale_fill_manual(values=colors)+
	theme_bw()+
	stat_compare_means(comparisons=Class,label="p.signif")+
	theme(axis.text.x=element_text(colour="black",size=14),	#设置x轴刻度标签的字体属性
		axis.text.y=element_text(size=12,face="plain"), 		#设置x轴刻度标签的字体属性
		axis.title.y=element_text(size = 15,face="plain"),	#设置y轴的标题的字体属性

		plot.title = element_text(size=12,face="bold",hjust = 0.5), #设置总标题的字体属性
		panel.grid.major = element_blank()
)

pdf('/public/workspace/yumiao/TE_GAST/CytoTRACE_score.pdf',width=5,height=5)
print(p1)
dev.off()
}

