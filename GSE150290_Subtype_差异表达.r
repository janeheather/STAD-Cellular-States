
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

set.seed(520)
STAD<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/STAD_P_E_T_N.rds')
DefaultAssay(STAD) <- "RNA"

STAD$MP_Subtype='Other'
STAD$MP_Subtype[which(STAD@meta.data$seurat_clusters %in% c(7,8))]='MP_1'
STAD$MP_Subtype[which(STAD@meta.data$seurat_clusters %in% c(5,9))]='MP_2'
STAD$MP_Subtype[which(STAD@meta.data$seurat_clusters %in% c(3,4,11))]='MP_3'

MP_Subtype=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")
pdf('/public/workspace/liuqzh/gastric_cancer/GSEA/GSE150290_Subtype_Dimplot.pdf',height=6,width=6)
DimPlot(STAD,group.by="MP_Subtype",reduction='tsne',pt.size=0.1,label = TRUE,cols=MP_Subtype)
DimPlot(STAD,group.by="MP_Subtype",reduction='umap',pt.size=0.1,label = TRUE,cols=MP_Subtype)
dev.off()

Idents(STAD)=STAD$MP_Subtype
STAD<-subset(STAD,subset=seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11))

{
marker=FindMarkers(object = STAD, ident.1=c("MP_1"),min.pct=0.3)
marker=marker[which(marker$p_val_adj < 0.05),]
####P值太显著,可以选取显著的基因，忽略掉这一值，用Percentage Difference来作为横坐标###

library(ggrepel)
###构建中间数据
dif=as.data.frame(marker)
dif$Difference=""
dif$Difference=(dif$pct.1-dif$pct.2)*100

###自己定义上下调基因
# (1) define up and down
log2FC=1
padj=0.01
dif$threshold="ns"
dif[which(dif$avg_log2FC > log2FC & dif$p_val_adj <padj & dif$Difference > 50 ),]$threshold="up"
dif[which(dif$avg_log2FC < (-log2FC) & dif$p_val_adj <padj & dif$Difference < -10 ),]$threshold="down"
dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
table(dif$threshold)

##选择可视化的前20个基因
up.genes <- rownames(dif[which(dif$threshold == "up")[1:10],])
down.genes <- rownames(dif[which(dif$threshold == "down")[1],])
diftop20gene<-c(as.character(up.genes),as.character(down.genes))
dif$lable=""
dif$lable[which(rownames(dif) %in% diftop20gene)]<-diftop20gene

pdf('/public/workspace/liuqzh/gastric_cancer/GSEA/GSE150290_Subtype_MP1_marker_gene.pdf',height=5,width=5)
ggplot(dif,aes(x=Difference,y=avg_log2FC))+
	geom_point(size=0.5,aes(color=threshold))+
	scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
	geom_label_repel(data=dif,aes(label=lable),label.padding = 0.1,fill="#F4F0EB",
	segment.size = 0.25,size=3,max.overlaps=100,force = 10)+
	geom_vline(xintercept = 0,linetype="dashed")+
	geom_hline(yintercept = 0,linetype="dashed")+
	labs(x="Percentage Difference",x="Avg_log2FC",title="MP1 VS Other")+
	theme_classic()
dev.off()
}

{
marker=FindMarkers(object = STAD, ident.1=c("MP_2"),min.pct=0.3)
marker=marker[which(marker$p_val_adj < 0.05),]
####P值太显著,可以选取显著的基因，忽略掉这一值，用Percentage Difference来作为横坐标###

library(ggrepel)
###构建中间数据
dif=as.data.frame(marker)
dif$Difference=""
dif$Difference=(dif$pct.1-dif$pct.2)*100

###自己定义上下调基因
# (1) define up and down
log2FC=1
padj=0.01
dif$threshold="ns"
dif[which(dif$avg_log2FC > log2FC & dif$p_val_adj <padj & dif$Difference > 40 ),]$threshold="up"
dif[which(dif$avg_log2FC < (-log2FC) & dif$p_val_adj <padj & dif$Difference < -10 ),]$threshold="down"
dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
table(dif$threshold)

##选择可视化的前20个基因
up.genes <- rownames(dif[which(dif$threshold == "up")[1:10],])
down.genes <- rownames(dif[which(dif$threshold == "down")[1:10],])
diftop20gene<-c(as.character(up.genes),as.character(down.genes))
dif$lable=""
dif$lable[which(rownames(dif) %in% diftop20gene)]<-diftop20gene

pdf('/public/workspace/liuqzh/gastric_cancer/GSEA/GSE150290_Subtype_MP2_marker_gene.pdf',height=5,width=5)
ggplot(dif,aes(x=Difference,y=avg_log2FC))+
	geom_point(size=0.5,aes(color=threshold))+
	scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
	geom_label_repel(data=dif,aes(label=lable),label.padding = 0.1,fill="#F4F0EB",
	segment.size = 0.25,size=3,max.overlaps=100,force = 10)+
	geom_vline(xintercept = 0,linetype="dashed")+
	geom_hline(yintercept = 0,linetype="dashed")+
	labs(x="Percentage Difference",x="Avg_log2FC",title="MP2 VS Other")+
	theme_classic()
dev.off()
}

{
marker=FindMarkers(object = STAD, ident.1=c("MP_3"),min.pct=0.3)
marker=marker[which(marker$p_val_adj < 0.05),]
####P值太显著,可以选取显著的基因，忽略掉这一值，用Percentage Difference来作为横坐标###

library(ggrepel)
###构建中间数据
dif=as.data.frame(marker)
dif$Difference=""
dif$Difference=(dif$pct.1-dif$pct.2)*100

###自己定义上下调基因
# (1) define up and down
log2FC=1
padj=0.01
dif$threshold="ns"
dif[which(dif$avg_log2FC > log2FC & dif$p_val_adj <padj & dif$Difference > 10 ),]$threshold="up"
dif[which(dif$avg_log2FC < (-log2FC) & dif$p_val_adj <padj & dif$Difference < -10 ),]$threshold="down"
dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
table(dif$threshold)

##选择可视化的前20个基因
up.genes <- rownames(dif[which(dif$threshold == "up")[1:10],])
down.genes <- rownames(dif[which(dif$threshold == "down")[1:10],])
diftop20gene<-c(as.character(up.genes),as.character(down.genes))
dif$lable=""
dif$lable[which(rownames(dif) %in% diftop20gene)]<-diftop20gene

pdf('/public/workspace/liuqzh/gastric_cancer/GSEA/GSE150290_Subtype_MP3_marker_gene.pdf',height=5,width=5)
ggplot(dif,aes(x=Difference,y=avg_log2FC))+
	geom_point(size=0.5,aes(color=threshold))+
	scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
	geom_label_repel(data=dif,aes(label=lable),label.padding = 0.1,fill="#F4F0EB",
	segment.size = 0.25,size=3,max.overlaps=100,force = 10)+
	geom_vline(xintercept = 0,linetype="dashed")+
	geom_hline(yintercept = 0,linetype="dashed")+
	labs(x="Percentage Difference",x="Avg_log2FC",title="MP3 VS Other")+
	theme_classic()
dev.off()
}

