
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
library(scrabble)
library(stringr)

set.seed(520)
STAD_P<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
STAD_MP<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')

STAD_P<-subset(STAD_P,subset=orig.ident %in% c('P01_EGC','P02_AGC','P04_EGC','P05_EGC','P06_EGC','P07_EGC',
											   'P08_EGC','P09_EGC','P10_EGC','P11_AGC','P12_AGC','P13_EGC',
											   'P15_AGC','P16_EGC','P17_EGC','P18_EGC','P19_EGC','P20_AGC',
											   'P21_EGC','P22_AGC','P23_AGC'))

STAD_MP<-subset(STAD_MP,subset=dataset %in% c('GSE150290'))

STAD_P@meta.data$cibresort_subtype=STAD_P@meta.data$celltype
STAD_P@meta.data$cibresort_subtype[which(rownames(STAD_P@meta.data) %in% rownames(STAD_MP@meta.data[which(STAD_MP@meta.data$Subtype == 'MP_1'),]))]='MP_1'
STAD_P@meta.data$cibresort_subtype[which(rownames(STAD_P@meta.data) %in% rownames(STAD_MP@meta.data[which(STAD_MP@meta.data$Subtype == 'MP_2'),]))]='MP_2'
STAD_P@meta.data$cibresort_subtype[which(rownames(STAD_P@meta.data) %in% rownames(STAD_MP@meta.data[which(STAD_MP@meta.data$Subtype == 'MP_3'),]))]='MP_3'
Idents(STAD_P)<-STAD_P$cibresort_subtype

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
subset_cells_type6 <- STAD_P[, Idents(STAD_P) %in% "Fibroblast"]
subset_cells_type6 <- subset_cells_type6[, sample(1:ncol(subset_cells_type6), size = round(0.1 * ncol(subset_cells_type6)))]
subset_cells_type7 <- STAD_P[, Idents(STAD_P) %in% "Endothelial"]
subset_cells_type7 <- subset_cells_type7[, sample(1:ncol(subset_cells_type7), size = round(0.1 * ncol(subset_cells_type7)))]
subset_cells_type8 <- STAD_P[, Idents(STAD_P) %in% "Endocrine"]
subset_cells_type8 <- subset_cells_type8[, sample(1:ncol(subset_cells_type8), size = round(0.1 * ncol(subset_cells_type8)))]

subset_cells_type9 <- STAD_P[, Idents(STAD_P) %in% "MP_1"]
subset_cells_type9 <- subset_cells_type9[, sample(1:ncol(subset_cells_type9), size = round(0.1 * ncol(subset_cells_type9)))]
subset_cells_type10 <- STAD_P[, Idents(STAD_P) %in% "MP_2"]
subset_cells_type10 <- subset_cells_type10[, sample(1:ncol(subset_cells_type10), size = round(0.1 * ncol(subset_cells_type10)))]
subset_cells_type11 <- STAD_P[, Idents(STAD_P) %in% "MP_3"]
subset_cells_type11 <- subset_cells_type11[, sample(1:ncol(subset_cells_type11), size = round(0.1 * ncol(subset_cells_type11)))]

###########
expr1<-as.data.frame(subset_cells_type1@assays$RNA@data)
colnames(expr1)=rep('B cell', ncol(expr1))

expr2<-as.data.frame(subset_cells_type2@assays$RNA@data)
colnames(expr2)=rep('T cell', ncol(expr2))

expr3<-as.data.frame(subset_cells_type3@assays$RNA@data)
colnames(expr3)=rep('NK cell', ncol(expr3))

expr4<-as.data.frame(subset_cells_type4@assays$RNA@data)
colnames(expr4)=rep('Macrophage', ncol(expr4))

expr5<-as.data.frame(subset_cells_type5@assays$RNA@data)
colnames(expr5)=rep('Mast cell', ncol(expr5))

expr6<-as.data.frame(subset_cells_type6@assays$RNA@data)
colnames(expr6)=rep('Fibroblast', ncol(expr6))

expr7<-as.data.frame(subset_cells_type7@assays$RNA@data)
colnames(expr7)=rep('Endothelial', ncol(expr7))

expr8<-as.data.frame(subset_cells_type8@assays$RNA@data)
colnames(expr8)=rep('Endocrine', ncol(expr8))

expr9<-as.data.frame(subset_cells_type9@assays$RNA@data)
colnames(expr9)=rep('MP_1', ncol(expr9))

expr10<-as.data.frame(subset_cells_type10@assays$RNA@data)
colnames(expr10)=rep('MP_2', ncol(expr10))

expr11<-as.data.frame(subset_cells_type11@assays$RNA@data)
colnames(expr11)=rep('MP_3', ncol(expr11))

sc_cibresort_ref=cbind(rownames(expr1),expr1,expr2,expr3,expr4,expr5,
					   expr6,expr7,expr8,expr9,expr10,
					   expr11)
colnames(sc_cibresort_ref)[1]='Gene'

write.table(sc_cibresort_ref,'/public/workspace/liuqzh/gastric_cancer/cibersort/GSE150290_sc_cibresort_ref.txt',sep='\t',row.names=F)

