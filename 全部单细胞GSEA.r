
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
a<-arrange(marker,desc(marker[,2]))

fc_MP_1<-a[,2]
names(fc_MP_1)<-rownames(a)

hallmark <- read.gmt("/public/workspace/yumiao/data/1209/h.all.v2022.1.Hs.symbols.gmt")

gsea.re1<- GSEA(fc_MP_1,  #待富集的基因列表
    TERM2GENE = hallmark,  #基因集
    pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'fdr',eps=0)  #指定 p 值校正方法
     
ss=as.data.frame(gsea.re1)
ss[,1:5]
GSEA_MP_1=ss[,1:8]
saveRDS(GSEA_MP_1,'/public/workspace/yumiao/STAD/gsea/fig_gsea/GSEA_MP_1.rds')

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_1_E2F_TARGETS.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'E2F_TARGETS',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.75,pvalY = 0.85,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_1_TNFA_SIGNALING_VIA_NFKB.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'TNFA_SIGNALING_VIA_NFKB',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.05,pvalY = 0.15,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_1_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf',height=4,width=7)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'TNFA_SIGNALING_VIA_NFKB',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.05,pvalY = 0.15,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()
}

############
############
{
marker=FindMarkers(object = STAD, ident.1=c("NES_MP_2"),min.pct=0.3)
a<-arrange(marker,desc(marker[,2]))

fc_MP_2<-a[,2]
names(fc_MP_2)<-rownames(a)

hallmark <- read.gmt("/public/workspace/yumiao/data/1209/h.all.v2022.1.Hs.symbols.gmt")

gsea.re1<- GSEA(fc_MP_2,  #待富集的基因列表
    TERM2GENE = hallmark,  #基因集
    pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'fdr',eps=0)  #指定 p 值校正方法
     
ss=as.data.frame(gsea.re1)
ss[,1:5]
GSEA_MP_2=ss[,1:8]
saveRDS(GSEA_MP_2,'/public/workspace/yumiao/STAD/gsea/fig_gsea/GSEA_MP_2.rds')

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_2_APICAL_JUNCTION.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'APICAL_JUNCTION',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.75,pvalY = 0.85,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_2_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'EPITHELIAL_MESENCHYMAL_TRANSITION',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.75,pvalY = 0.85,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_2_TNFA_SIGNALING_VIA_NFKB.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'TNFA_SIGNALING_VIA_NFKB',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.05,pvalY = 0.15,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_2_E2F_TARGETS.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'E2F_TARGETS',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.05,pvalY = 0.15,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()
}

############
############
{
marker=FindMarkers(object = STAD, ident.1=c("NES_MP_3"),min.pct=0)
a<-arrange(marker,desc(marker[,2]))

fc_MP_3<-a[,2]
names(fc_MP_3)<-rownames(a)

hallmark <- read.gmt("/public/workspace/yumiao/data/1209/h.all.v2022.1.Hs.symbols.gmt")

gsea.re1<- GSEA(fc_MP_3,  #待富集的基因列表
    TERM2GENE = hallmark,  #基因集
    pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'fdr',eps=0)  #指定 p 值校正方法
     
ss=as.data.frame(gsea.re1)
GSEA_MP_3=ss[,1:8]
saveRDS(GSEA_MP_3,'/public/workspace/yumiao/STAD/gsea/fig_gsea/GSEA_MP_3.rds')

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_3_TNFA_SIGNALING_VIA_NFKB.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'TNFA_SIGNALING_VIA_NFKB',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.75,pvalY = 0.85,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_3_HYPOXIA.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'HYPOXIA',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.75,pvalY = 0.85,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_3_INFLAMMATORY_RESPONSE.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'INFLAMMATORY_RESPONSE',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.75,pvalY = 0.85,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_3_OXIDATIVE_PHOSPHORYLATION.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'OXIDATIVE_PHOSPHORYLATION',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.05,pvalY = 0.15,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()

pdf('/public/workspace/yumiao/STAD/gsea/fig_gsea/MP_3_APICAL_JUNCTION.pdf',height=4,width=8)
p = gseaNb(object = gsea.re1,
		   geneSetID = 'APICAL_JUNCTION',
		   newGsea = T,addPoint = F,addPval = T,
		   pvalX = 0.05,pvalY = 0.15,
		   pCol = 'black',
		   pHjust = 0,rmPrefix = F,newCurveCol=c("#D1832C","#D1832C","#D1832C"),newHtCol=c("#BD3934","white","#5B93D4")
		   )
print(p)
dev.off()
}


############
library(ComplexHeatmap)
library(circlize)

GSEA_MP_1<-readRDS('/public/workspace/yumiao/STAD/gsea/fig_gsea/GSEA_MP_1.rds')
GSEA_MP_2<-readRDS('/public/workspace/yumiao/STAD/gsea/fig_gsea/GSEA_MP_2.rds')
GSEA_MP_3<-readRDS('/public/workspace/yumiao/STAD/gsea/fig_gsea/GSEA_MP_3.rds')

GSEA_MP_1=GSEA_MP_1[,5:6]
GSEA_MP_2=GSEA_MP_2[,5:6]
GSEA_MP_3=GSEA_MP_3[,5:6]

colnames(GSEA_MP_1)=c('MP_1_NES','MP_1_pvalue')
colnames(GSEA_MP_2)=c('MP_2_NES','MP_2_pvalue')
colnames(GSEA_MP_3)=c('MP_3_NES','MP_3_pvalue')

GSEA_MP_1=GSEA_MP_1[intersect(rownames(GSEA_MP_1),rownames(GSEA_MP_2)),]
GSEA_MP_2=GSEA_MP_2[intersect(rownames(GSEA_MP_1),rownames(GSEA_MP_2)),]
aa=cbind(GSEA_MP_1,GSEA_MP_2)

aa=aa[intersect(rownames(aa),rownames(GSEA_MP_3)),]
GSEA_MP_3=GSEA_MP_3[intersect(rownames(aa),rownames(GSEA_MP_3)),]

GSEA_plot_data=cbind(aa,GSEA_MP_3)

# 分离NES和pvalue
data=GSEA_plot_data

# 分离NES和pvalue，并转置
nes_data <- t(data[c("MP_1_NES", "MP_2_NES", "MP_3_NES")])  # 以此类推...
pvalue_data <- t(data[c("MP_1_pvalue", "MP_2_pvalue", "MP_3_pvalue")])  # 以此类推...

# 定义一个函数来添加星号
add_stars <- function(j, i, x, y, width, height) {
  if (pvalue_data[i, j] < 0.05) {
    pushViewport(viewport(x = x, y = y, width = width, height = height))
    grid.text("*", gp = gpar(col = "black", fontsize = 10))
    popViewport()
  }
}

# 绘制热图
ht <- Heatmap(nes_data,
              name = "NES",
              column_names_side = "bottom",
              show_row_names = TRUE,
              show_column_names = TRUE,
              clustering_distance_rows = "euclidean",
              clustering_distance_columns = "euclidean",
              clustering_method_rows = "complete",
              clustering_method_columns = "complete",
              col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              top_annotation = HeatmapAnnotation(pvalue = anno_mark(at = c(0.05), labels = "*", labels_gp = gpar(fontsize = 10, col = "red")))
)
ht
