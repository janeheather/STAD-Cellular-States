
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(dplyr)

STAD_Harmony<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')

STAD_Harmony$dataset<-'unrecognized'
STAD_Harmony$dataset[which(STAD_Harmony$orig.ident %in% c('t1_5846','t1_5866','t1_5931','t1_6207','t1_6342','t1_6592','t1_6709','t2_5866','t2_5931'))]<-'Anuja_data'
STAD_Harmony$dataset[which(STAD_Harmony$orig.ident %in% c('CRA002586-D01-T','CRA002586-D02-T','CRA002586-D03-T','CRA002586-P01-T',
														  'CRA002586-P02-T','CRA002586-P03-T','CRA002586-P04-T','CRA002586-P05-T',
														  'CRA002586-P06-T'))]<-'CRA002586'
STAD_Harmony$dataset[which(STAD_Harmony$orig.ident %in% c('P01_EGC','P02_AGC','P04_EGC','P05_EGC','P06_EGC','P07_EGC',
														  'P08_EGC','P09_EGC','P10_EGC','P11_AGC','P12_AGC','P13_EGC',
														  'P15_AGC','P16_EGC','P17_EGC','P18_EGC','P19_EGC','P20_AGC',
														  'P21_EGC','P22_AGC','P23_AGC'))]<-'GSE150290'
STAD_Harmony$dataset[which(STAD_Harmony$orig.ident %in% c('GSE167297-P1-1','GSE167297-P1-2','GSE167297-P2-1','GSE167297-P2-2',
														  'GSE167297-P3-1','GSE167297-P3-2','GSE167297-P4-1','GSE167297-P4-2'))]<-'GSE167297'
STAD_Harmony$dataset[which(STAD_Harmony$orig.ident %in% c('Tu1','Tu10','Tu11','Tu12','Tu13',
														  'Tu14','Tu15','Tu16','Tu17','Tu18',
														  'Tu19','Tu2','Tu20','Tu21','Tu22',
														  'Tu23','Tu24','Tu25','Tu26','Tu3',
														  'Tu4','Tu5','Tu6','Tu7','Tu8','Tu9'))]<-'GSE183904'

STAD_Harmony<-subset(STAD_Harmony,subset=celltype %in% c('Epithelium'))

list <- SplitObject(STAD_Harmony, split.by = "dataset")
for(i in 1:length(list)){
	tmp <- NormalizeData(list[[i]])
	tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
	list[[i]] <- tmp
}
anchors <- FindIntegrationAnchors(object.list = list, dims = 1:20)
STAD_Harmony <- IntegrateData(anchorset = anchors, dims = 1:20)
STAD_Harmony <- ScaleData(STAD_Harmony, features = rownames(STAD_Harmony))
STAD_Harmony <- RunPCA(STAD_Harmony, features = VariableFeatures(object = STAD_Harmony))
STAD_Harmony <- FindNeighbors(STAD_Harmony, dims = 1:20)
STAD_Harmony <- FindClusters(STAD_Harmony, resolution = 1)
STAD_Harmony <- RunUMAP(STAD_Harmony, dims = 1:20)
STAD_Harmony <- RunTSNE(STAD_Harmony, dims = 1:20)

DimPlot(STAD_Harmony,group.by="celltype",reduction='umap',label=F,pt.size=0.2)

#saveRDS(STAD_Harmony, "/public/workspace/liuqzh/gastric/STAD_Epithelium_Integrate_celltype.rds")

dat<-readRDS('/public/workspace/liuqzh/gastric/STAD_Epithelium_Integrate_celltype.rds')
DefaultAssay(dat)<-'RNA'

pdf("/public/workspace/liuqzh/gastric_cancer/figure1D_ne.pdf")
col1=c('0'='#BB3732','15'='#C95F5B','3'='#D68784','5'='#E4AFAD',
	'17'='#A3855B','4'='#B59D7C','12'='#C8B69D','6'='#DACEBD',
	'18'='#5991D2','7'='#7AA7DB','1'='#9BBDE4','14'='#BDD3ED',
	'21'='#CF812A','13'='#D99A55','22'='#E2B37F','20'='#ECCDAA',
	'11'='#509286','23'='#73A89E','16'='#96BEB6','30'='#B9D3CF',
	'9'='#427F39','8'='#689961','10'='#8EB288','28'='#B3CCB0',
	'29'='#282278','19'='#534E93','24'='#7E7AAE','31'='#7E7AAE',
	'2'='#9F71A7','25'='#B28DB9','26'='#C5AACA','27'='#D9C6DC','28'='#7E7AAE')
DimPlot(dat, reduction = "tsne",label=F,cols=col1)
dev.off()

DimPlot(dat, reduction = "tsne",label=T,cols=col1)

pdf("/public/workspace/liuqzh/gastric_cancer/figureS1A_ne.pdf")
col2=c('Anuja_data'='#BB3732',
	'CRA002586'='#A3855B',
	'GSE150290'='#5991D2',
	'GSE167297'='#CF812A',
	'GSE183904'='#509286')
DimPlot(dat, reduction = "tsne",group.by='dataset',label=F,cols=col2)
dev.off()

pdf("/public/workspace/liuqzh/gastric_cancer/figureS1F_ne.pdf",width=20,height=8)
DimPlot(dat, reduction = "tsne",group.by='orig.ident',label=F)
dev.off()

##############
bytlib load R-4.0.2
bytlib load gcc
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
library(dplyr)

STAD_P<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(STAD_P)<-'RNA'
Idents(STAD_P)=STAD_P$celltype
STAD_P=ScaleData(STAD_P,features = rownames(STAD_P))

Celltype.markers=FindAllMarkers(object = STAD_P, only.pos = TRUE, logfc.threshold = 0.3)
#saveRDS(Celltype.markers,"/public/workspace/liuqzh/gastric_cancer/Celltype.markers_gene_marker.rds")

Celltype.markers<-readRDS('/public/workspace/liuqzh/gastric_cancer/Celltype.markers_gene_marker.rds')
top10 <- Celltype.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)

colors=c("Epithelium"="#D1832C","Fibroblast"="#529387","Endothelial"="#5B93D4",
		 "NK cell"="#2A247A","T cell"="#A173A9","B cell"="#BD3833",
		 "Macrophage"="#751C51","Mast cell"="#44813B","Endocrine cell"="#A3875B")

# Downsample the number of cells per identity class
aaa <- subset(x = STAD_P, downsample = 1000)

aaa$celltype=factor(x=aaa$celltype,levels=c("B cell","T cell","NK cell","Macrophage","Mast cell","Epithelium","Fibroblast","Endothelial","Endocrine cell"))
DefaultAssay(aaa)<-'RNA'
pdf("/public/workspace/liuqzh/gastric_cancer/DoHeatmap.pdf",width=8,height=8)
DoHeatmap(aaa,   
          features = as.character(unique(top10$gene)),   
          group.by = "celltype",  
		  group.colors=colors)+
		  scale_fill_gradientn(colors = c("navy","white","firebrick3"))
dev.off()
###########
###########
library(ComplexHeatmap)
library(Seurat)
##提取标准化表达矩阵
#⚠️提取scale.data矩阵的时候一定要注意，做ScaleData()的时候一定是scale了所有的基因，而不是默认的2000个基因 
Celltype.markers<-readRDS('/public/workspace/liuqzh/gastric_cancer/Celltype.markers_gene_marker.rds')
STAD_P<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')

Idents(STAD_P)=STAD_P$celltype
#aaa <- subset(x = STAD_P, downsample = 1000)
table(STAD_P$celltype)

top_n=40
B_cell_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='B cell'),])[1:top_n]
T_cell_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='T cell'),])[1:top_n]
NK_cell_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='NK cell'),])[1:top_n]
Macrophage_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='Macrophage'),])[1:top_n]
Mast_cell_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='Mast cell'),])[1:top_n]
Epithelium_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='Epithelium'),])[1:top_n]
Fibroblast_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='Fibroblast'),])[1:top_n]
Endothelial_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='Endothelial'),])[1:top_n]
Endocrine_marker=rownames(Celltype.markers[which(Celltype.markers$cluster=='Endocrine'),])[1:top_n]
top10_gene=c(B_cell_marker,T_cell_marker,NK_cell_marker,
			Macrophage_marker,Mast_cell_marker,Epithelium_marker,
			Fibroblast_marker,Endothelial_marker,Endocrine_marker)

# 从细胞类型中分别取出10%的细胞
subset_cells_type1 <- STAD_P[, Idents(STAD_P) %in% "B cell"]
subset_cells_type1 <- subset_cells_type1[, sample(1:ncol(subset_cells_type1), size = round(0.1 * ncol(subset_cells_type1)))]
subset_cells_type2 <- STAD_P[, Idents(STAD_P) %in% "T cell"]
subset_cells_type2 <- subset_cells_type2[, sample(1:ncol(subset_cells_type2), size = round(0.1 * ncol(subset_cells_type2)))]
subset_cells_type3 <- STAD_P[, Idents(STAD_P) %in% "NK cell"]
subset_cells_type3 <- subset_cells_type3[, sample(1:ncol(subset_cells_type3), size = round(0.1 * ncol(subset_cells_type3)))]
subset_cells_type4 <- STAD_P[, Idents(STAD_P) %in% "Macrophage"]
subset_cells_type4 <- subset_cells_type4[, sample(1:ncol(subset_cells_type4), size = round(0.1 * ncol(subset_cells_type4)))]
subset_cells_type5 <- STAD_P[, Idents(STAD_P) %in% "Mast cell"]
subset_cells_type5 <- subset_cells_type5[, sample(1:ncol(subset_cells_type5), size = round(0.1 * ncol(subset_cells_type5)))]
subset_cells_type6 <- STAD_P[, Idents(STAD_P) %in% "Epithelium"]
subset_cells_type6 <- subset_cells_type6[, sample(1:ncol(subset_cells_type6), size = round(0.1 * ncol(subset_cells_type6)))]
subset_cells_type7 <- STAD_P[, Idents(STAD_P) %in% "Fibroblast"]
subset_cells_type7 <- subset_cells_type7[, sample(1:ncol(subset_cells_type7), size = round(0.1 * ncol(subset_cells_type7)))]
subset_cells_type8 <- STAD_P[, Idents(STAD_P) %in% "Endothelial"]
subset_cells_type8 <- subset_cells_type8[, sample(1:ncol(subset_cells_type8), size = round(0.1 * ncol(subset_cells_type8)))]
subset_cells_type9 <- STAD_P[, Idents(STAD_P) %in% "Endocrine"]
subset_cells_type9 <- subset_cells_type9[, sample(1:ncol(subset_cells_type9), size = round(0.1 * ncol(subset_cells_type9)))]

# 添加标识符列
subset_cells_type1$group <- "B cell"
subset_cells_type2$group <- "T cell"
subset_cells_type3$group <- "NK cell"
subset_cells_type4$group <- "Macrophage"
subset_cells_type5$group <- "Mast cell"
subset_cells_type6$group <- "Epithelium"
subset_cells_type7$group <- "Fibroblast"
subset_cells_type8$group <- "Endothelial"
subset_cells_type9$group <- "Endocrine"
# 合并子集
subset_data<-merge(subset_cells_type1,y= c(subset_cells_type2,subset_cells_type3,
				   subset_cells_type4,subset_cells_type5,subset_cells_type6,
				   subset_cells_type7,subset_cells_type8,subset_cells_type9))
subset_data=ScaleData(subset_data,features = rownames(subset_data))
#saveRDS(subset_data,"/public/workspace/liuqzh/gastric_cancer/subset_data_ScaleData.rds")

mat <- GetAssayData(subset_data,slot = 'scale.data')
##获得基因和细胞聚类信息
gene_features <- top10_gene
cluster_info <- subset_data$group
##筛选矩阵
mat <- as.matrix(mat[intersect(gene_features,rownames(mat)),names(cluster_info)])
##输入想在图上展示出来的marker基因，获得基因在热图中的位置信息
gene <- c('CD3E','GNLY','GZMB','IGJ','IGLL5','JCHAIN')
gene_pos <- which(rownames(mat)%in%gene)
row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,labels = gene))
##画个热图看看
Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        right_annotation = row_anno)

col <- c("Epithelium"="#D1832C","Fibroblast"="#529387","Endothelial"="#5B93D4",
		 "NK cell"="#2A247A","T cell"="#A173A9","B cell"="#BD3833",
		 "Macrophage"="#751C51","Mast cell"="#44813B","Endocrine"="#A3875B")
names(col) <- levels(cluster_info)

top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=col),
                                                 labels = levels(cluster_info),
                                                 labels_gp = gpar(cex=0.5,col='white')))
library(circlize)
col_fun = colorRamp2(c(-2, 1, 4), c("#377EB8", "white", "#E41A1C"))

Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
		show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno, #在热图边上增加注释	
		column_title = NULL,
        right_annotation = row_anno,
        heatmap_legend_param = list(
          title='Expression',
          title_position='leftcenter-rot'),
        col = col_fun)
		
#######################
bytlib load R-4.0.2
bytlib load gcc
R

library(Seurat)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(scran)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)

Celltype.markers<-readRDS('/public/workspace/liuqzh/gastric_cancer/Celltype.markers_gene_marker.rds')
STAD_P<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
Celltype.markers=Celltype.markers[which(Celltype.markers$pct.1>0.5),]
top5 <- Celltype.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
Idents(STAD_P)<-STAD_P$celltype

p1 <- DotPlot(STAD_P,
    features = split(top5$gene, top5$cluster),
    cols = c("grey", '#BD3833')
    ) +
    RotatedAxis() + # 来自Seurat
    theme(
         panel.border = element_rect(color = "black"),
         panel.spacing = unit(1, "mm"),
         axis.title = element_blank(),
         axis.text.y = element_blank(),
     )
     p1

pdf("/public/workspace/liuqzh/gastric_cancer/Dotplot.pdf",width=13,height=4)
p1
dev.off()
