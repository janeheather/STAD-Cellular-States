#########PLF_EMT_NLT的反卷积###########
bytlib load R-4.0.2
bytlib load gcc
R

library(e1071)
library(preprocessCore)
library(parallel)
library(CIBERSORT)
library(Seurat)
library(dplyr)

STAD_FPKM=readRDS("/public/workspace/liuqzh/gastric/STAD_FPKM_data_exp.rds")
STAD_FPKM=STAD_FPKM[,grep('01A',colnames(STAD_FPKM))]
colnames(STAD_FPKM)=gsub('-01A','',colnames(STAD_FPKM))
#write.table(STAD_FPKM, "/public/workspace/liuqzh/gastric_cancer/cibersort/STAD_FPKM_cibersort.txt", row.names = T, sep = "\t")

genelist<-read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)

STAD_Integrate_celltype=readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
STAD_MP_Subtype<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')

STAD_Integrate_celltype$cibersort_subtype=STAD_Integrate_celltype$celltype
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_MP_Subtype@meta.data[which(STAD_MP_Subtype@meta.data$Subtype=='MP_1'),]))]<-'MP1'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_MP_Subtype@meta.data[which(STAD_MP_Subtype@meta.data$Subtype=='MP_2'),]))]<-'MP2'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_MP_Subtype@meta.data[which(STAD_MP_Subtype@meta.data$Subtype=='MP_3'),]))]<-'MP3'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_Integrate_celltype@meta.data[which(STAD_Integrate_celltype@meta.data$celltype=='Endothelial'),]))]<-'Stormal'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_Integrate_celltype@meta.data[which(STAD_Integrate_celltype@meta.data$celltype=='Fibroblast'),]))]<-'Stormal'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_Integrate_celltype@meta.data[which(STAD_Integrate_celltype@meta.data$celltype=='B cell'),]))]<-'Immune'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_Integrate_celltype@meta.data[which(STAD_Integrate_celltype@meta.data$celltype=='NK cell'),]))]<-'Immune'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_Integrate_celltype@meta.data[which(STAD_Integrate_celltype@meta.data$celltype=='T cell'),]))]<-'Immune'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_Integrate_celltype@meta.data[which(STAD_Integrate_celltype@meta.data$celltype=='Macrophage'),]))]<-'Immune'
STAD_Integrate_celltype$cibersort_subtype[which(colnames(STAD_Integrate_celltype) %in% rownames(STAD_Integrate_celltype@meta.data[which(STAD_Integrate_celltype@meta.data$celltype=='Mast cell'),]))]<-'Immune'

#DimPlot(STAD_Integrate_celltype,group.by="cibersort_subtype",reduction = "umap",pt.size=0.1,label=TRUE)
Idents(STAD_Integrate_celltype)<-STAD_Integrate_celltype$cibersort_subtype

averger_exp <- AverageExpression(STAD_Integrate_celltype)[[1]]
#markers_MP_1 <- FindMarkers(STAD_Integrate_celltype, ident.1 = "MP1",min.pct = 0.3)
#markers_MP_2 <- FindMarkers(STAD_Integrate_celltype, ident.1 = "MP2",min.pct = 0.3)
#markers_MP_3 <- FindMarkers(STAD_Integrate_celltype, ident.1 = "MP3",min.pct = 0.3)
#saveRDS(markers_MP_1,"/public/workspace/liuqzh/gastric_cancer/markers_MP_1.rds")
#saveRDS(markers_MP_2,"/public/workspace/liuqzh/gastric_cancer/markers_MP_2.rds")
#saveRDS(markers_MP_3,"/public/workspace/liuqzh/gastric_cancer/markers_MP_3.rds")

markers_MP_1=readRDS('/public/workspace/liuqzh/gastric_cancer/markers_MP_1.rds')
markers_MP_2=readRDS('/public/workspace/liuqzh/gastric_cancer/markers_MP_2.rds')
markers_MP_3=readRDS('/public/workspace/liuqzh/gastric_cancer/markers_MP_3.rds')

markers_MP_1=markers_MP_1[which(markers_MP_1$avg_log2FC > 0.5 & markers_MP_1$pct.2 < 0.5 & markers_MP_1$pct.1/markers_MP_1$pct.2 > 2),]
markers_MP_1=markers_MP_1[order(markers_MP_1$avg_log2FC,decreasing = TRUE),]
MP_1=rownames(markers_MP_1)[1:30]
markers_MP_2=markers_MP_2[which(markers_MP_2$avg_log2FC > 0.5 & markers_MP_2$pct.2 < 0.5 & markers_MP_2$pct.1/markers_MP_2$pct.2 > 2),]
markers_MP_2=markers_MP_2[order(markers_MP_2$avg_log2FC,decreasing = TRUE),]
MP_2=rownames(markers_MP_2)[1:30]
markers_MP_3=markers_MP_3[which(markers_MP_3$avg_log2FC > 0.5 & markers_MP_3$pct.2 < 0.5 & markers_MP_3$pct.1/markers_MP_3$pct.2 > 2),]
markers_MP_3=markers_MP_3[order(markers_MP_2$avg_log2FC,decreasing = TRUE),]
MP_3=rownames(markers_MP_3)[1:30]

Celltype_markers_gene<-readRDS('/public/workspace/liuqzh/gastric_cancer/Celltype.markers_gene_marker.rds')
gene_list=Celltype_markers_gene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% as.data.frame()

genes=c(gene_list$gene,MP_1,MP_2,MP_3) 
averger=averger_exp[which(rownames(averger_exp) %in% genes),]
sig_matrix=averger

#write.table(sig_matrix, "/public/workspace/liuqzh/gastric_cancer/cibersort/sig_matrix_top20_Stor_Immu.txt", row.names = T, sep = "\t")


rowscale <- results[,1:ncol(sig_matrix)]#只是相当于备份了一下results
rowscale <- rowscale[,apply(rowscale, 2, function(x){sum(x)>0})]#删除全是0的列
pheatmap(rowscale,
         scale = 'row',#按行标准化，不标准化就会按绝对值显示，很诡异
         cluster_col=T,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
         cluster_row=T,#是否对行聚类
         angle_col = "315")#调整X轴坐标的倾斜角度