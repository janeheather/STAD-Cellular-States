
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(monocle)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)
library(ggrepel)
library(psych)
library(pheatmap)
library(GseaVis)
library(dplyr)
library(clusterProfiler)
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library('GSEABase')
library(fgsea)
library(org.Hs.eg.db)
library(presto)

STAD<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(STAD) <- "RNA"

STAD@meta.data$NES_MP_Subtype='Intermediate'
STAD@meta.data$NES_MP_Subtype[which(STAD@meta.data$MP_1>0 & STAD@meta.data$Subtype == 'MP_1')]='NES_MP_1'
STAD@meta.data$NES_MP_Subtype[which(STAD@meta.data$MP_2>0 & STAD@meta.data$Subtype == 'MP_2')]='NES_MP_2'
STAD@meta.data$NES_MP_Subtype[which(STAD@meta.data$MP_3>0 & STAD@meta.data$Subtype == 'MP_3')]='NES_MP_3'

NES_MP_Subtype=c("NES_MP_1"="#D1832C","NES_MP_2"="#5B93D4","NES_MP_3"="#A97598")
pdf('/public/workspace/liuqzh/gastric_cancer/GSEA/Subtype_Dimplot.pdf',height=6,width=6)
DimPlot(STAD,group.by="NES_MP_Subtype",reduction='tsne',pt.size=0.1,label = TRUE,cols=NES_MP_Subtype)
DimPlot(STAD,group.by="NES_MP_Subtype",reduction='umap',pt.size=0.1,label = TRUE,cols=NES_MP_Subtype)
dev.off()

Idents(STAD)=STAD$NES_MP_Subtype

{
marker=FindMarkers(object = STAD, ident.1=c("NES_MP_1"),min.pct=0.3)
marker=marker[which(marker$p_val_adj < 0.05),]
####P值太显著,可以选取显著的基因，忽略掉这一值，用Percentage Difference来作为横坐标###

library(ggrepel)
###构建中间数据
dif=as.data.frame(marker)
dif$Difference=""
dif$Difference=(dif$pct.1-dif$pct.2)*100

###自己定义上下调基因
# (1) define up and down
log2FC=0.5
padj=0.01
dif$threshold="ns"
dif[which(dif$avg_log2FC > log2FC & dif$p_val_adj <padj & dif$Difference > 20 ),]$threshold="up"
dif[which(dif$avg_log2FC < (-log2FC) & dif$p_val_adj <padj & dif$Difference < -20 ),]$threshold="down"
dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
table(dif$threshold)

##选择可视化的前20个基因
up.genes <- head(dif$X[which(dif$threshold == "up")],20)
down.genes <- head(dif$X[which(dif$threshold=="down")],3)
dif$lable=""
diftop20gene<-c(as.character(up.genes),as.character(down.genes))
dif$lable[match(diftop20gene,dif$X)]<-diftop20gene

pdf('/public/workspace/liuqzh/gastric_cancer/GSEA/Subtype_MP1_marker_gene.pdf',height=5,width=5)
ggplot(dif,aes(x=Difference,y=avg_log2FC))+
	geom_point(size=0.5,aes(color=threshold))+
	scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
	geom_label_repel(data=dif,aes(label=lable),label.padding = 0.1,fill="tomato2", 
	segment.size = 0.25,size=5,max.overlaps=100,force = 10)+
	geom_vline(xintercept = 0,linetype="dashed")+
	geom_hline(yintercept = 0,linetype="dashed")+
	labs(x="Percentage Difference",x="Avg_log2FC",title="MP1 VS MP2/MP3")+
	theme_classic()
dev.off()	

}

{
marker=FindMarkers(object = STAD, ident.1=c("NES_MP_2"),min.pct=0.3)
marker=marker[which(marker$p_val_adj < 0.05),]
####P值太显著,可以选取显著的基因，忽略掉这一值，用Percentage Difference来作为横坐标###

library(ggrepel)
###构建中间数据
dif=as.data.frame(marker)
dif$Difference=""
dif$Difference=(dif$pct.1-dif$pct.2)*100

###自己定义上下调基因
# (1) define up and down
log2FC=0.5
padj=0.01
dif$threshold="ns"
dif[which(dif$avg_log2FC > log2FC & dif$p_val_adj <padj & dif$Difference > 20 ),]$threshold="up"
dif[which(dif$avg_log2FC < (-log2FC) & dif$p_val_adj <padj & dif$Difference < -20 ),]$threshold="down"
dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
table(dif$threshold)

##选择可视化的前20个基因
up.genes <- head(dif$X[which(dif$threshold == "up")],20)
down.genes <- head(dif$X[which(dif$threshold=="down")],20)
dif$lable=""
diftop20gene<-c(as.character(up.genes),as.character(down.genes))
dif$lable[match(diftop20gene,dif$X)]<-diftop20gene

pdf('/public/workspace/liuqzh/gastric_cancer/GSEA/Subtype_MP2_marker_gene.pdf',height=5,width=5)
ggplot(dif,aes(x=Difference,y=avg_log2FC))+
	geom_point(size=0.5,aes(color=threshold))+
	scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
	geom_label_repel(data=dif,aes(label=lable),label.padding = 0.1,fill="tomato2", 
	segment.size = 0.25,size=5,max.overlaps=100,force = 10)+
	geom_vline(xintercept = 0,linetype="dashed")+
	geom_hline(yintercept = 0,linetype="dashed")+
	labs(x="Percentage Difference",x="Avg_log2FC",title="MP2 VS MP1/MP3")+
	theme_classic()
dev.off()	
}

{
marker=FindMarkers(object = STAD, ident.1=c("NES_MP_3"),min.pct=0.3)
marker=marker[which(marker$p_val_adj < 0.05),]
####P值太显著,可以选取显著的基因，忽略掉这一值，用Percentage Difference来作为横坐标###

library(ggrepel)
###构建中间数据
dif=as.data.frame(marker)
dif$Difference=""
dif$Difference=(dif$pct.1-dif$pct.2)*100

###自己定义上下调基因
# (1) define up and down
log2FC=0.5
padj=0.01
dif$threshold="ns"
dif[which(dif$avg_log2FC > log2FC & dif$p_val_adj <padj & dif$Difference > 20 ),]$threshold="up"
dif[which(dif$avg_log2FC < (-log2FC) & dif$p_val_adj <padj & dif$Difference < -20 ),]$threshold="down"
dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
table(dif$threshold)

##选择可视化的前20个基因
up.genes <- head(dif$X[which(dif$threshold == "up")],15)
down.genes <- head(dif$X[which(dif$threshold=="down")],20)
dif$lable=""
diftop20gene<-c(as.character(up.genes),as.character(down.genes))
dif$lable[match(diftop20gene,dif$X)]<-diftop20gene

pdf('/public/workspace/liuqzh/gastric_cancer/GSEA/Subtype_MP3_marker_gene.pdf',height=5,width=5)
ggplot(dif,aes(x=Difference,y=avg_log2FC))+
	geom_point(size=0.5,aes(color=threshold))+
	scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
	geom_label_repel(data=dif,aes(label=lable),label.padding = 0.1,fill="tomato2", 
	segment.size = 0.25,size=5,max.overlaps=100,force = 10)+
	geom_vline(xintercept = 0,linetype="dashed")+
	geom_hline(yintercept = 0,linetype="dashed")+
	labs(x="Percentage Difference",x="Avg_log2FC",title="MP3 VS MP1/MP2")+
	theme_classic()
dev.off()
}









