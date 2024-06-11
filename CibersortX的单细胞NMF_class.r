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

genelist<-read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv",sep = ",",header = T)

STAD_Integrate_celltype=readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
STAD_MP_Subtype<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')

scdata_score<-readRDS('/public/workspace/liuqzh/gastric_cancer/scalop-master/scdata_score.rds')
Com1=c('Mac_MP_3','Endo_MP_1','Endo_MP_7','Endo_MP_5','Fib_MP_1',
	   'Fib_MP_4','Maligant_MP_2','Mast_MP_4','B_MP_6')
Com2=c('B_MP_5','T_cell_MP_1','Mast_MP_3','Endo_MP_6','Maligant_MP_1',
	   'T_cell_MP_3','B_MP_4','T_cell_MP_4','B_MP_3','Mac_MP_2','Endo_MP_3')
Com3=c('Mac_MP_1','Fib_MP_2')
Com4=c('Endo_MP_2','Mast_MP_2','Fib_MP_3','T_cell_MP_5','Maligant_MP_3',
	   'Mast_MP_1','T_cell_MP_6','B_MP_2','T_cell_MP_2','B_MP_1','Endo_MP_4')

STAD_Integrate_celltype@meta.data=STAD_Integrate_celltype@meta.data[intersect(rownames(scdata_score),rownames(STAD_Integrate_celltype@meta.data)),]
scdata_score=cbind(scdata_score,STAD_Integrate_celltype@meta.data$celltype)
colnames(scdata_score)[ncol(scdata_score)]='Celltype'
scdata_score=as.data.frame(scdata_score)

Fib=scdata_score[which(scdata_score$Celltype %in% 'Fibroblast'),grep('Fib_',colnames(scdata_score))]
Fib$NMF_Class=paste('Fib_MP_',apply(Fib,1,which.max),sep='')
Mac=scdata_score[which(scdata_score$Celltype %in% 'Macrophage'),grep('Mac_',colnames(scdata_score))]
Mac$NMF_Class=paste('Mac_MP_',apply(Mac,1,which.max),sep='')
T_cell=scdata_score[which(scdata_score$Celltype %in% c('NK cell','T cell')),grep('T_cell_',colnames(scdata_score))]
T_cell$NMF_Class=paste('T_cell_MP_',apply(T_cell,1,which.max),sep='')
Mast=scdata_score[which(scdata_score$Celltype %in% 'Mast cell'),grep('Mast_',colnames(scdata_score))]
Mast$NMF_Class=paste('Mast_MP_',apply(Mast,1,which.max),sep='')
B=scdata_score[which(scdata_score$Celltype %in% 'B cell'),grep('B_',colnames(scdata_score))]
B$NMF_Class=paste('B_MP_',apply(B,1,which.max),sep='')
Endo=scdata_score[which(scdata_score$Celltype %in% 'Endothelial'),grep('Endo_',colnames(scdata_score))]
Endo$NMF_Class=paste('Endo_MP_',apply(Endo,1,which.max),sep='')
Maligant=scdata_score[which(scdata_score$Celltype %in% 'Epithelium'),grep('Maligant_',colnames(scdata_score))]
Maligant$NMF_Class=paste('Maligant_MP_',apply(Maligant,1,which.max),sep='')

DATA = rbind(Fib[,ncol(Fib),drop=F],Mac[,ncol(Mac),drop=F],T_cell[,ncol(T_cell),drop=F],
			   Mast[,ncol(Mast),drop=F],B[,ncol(B),drop=F],Endo[,ncol(Endo),drop=F],
			   Maligant[,ncol(Maligant),drop=F])
DATA=as.data.frame(DATA)

STAD_Integrate_celltype=readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
STAD_Integrate_celltype$NMF_Class<-'unrecognized'
STAD_Integrate_celltype$NMF_Class[which(rownames(STAD_Integrate_celltype@meta.data) %in% rownames(DATA[which(DATA$NMF_Class %in% Com1),,drop=F]))]<-'Com1'
STAD_Integrate_celltype$NMF_Class[which(rownames(STAD_Integrate_celltype@meta.data) %in% rownames(DATA[which(DATA$NMF_Class %in% Com2),,drop=F]))]<-'Com2'
STAD_Integrate_celltype$NMF_Class[which(rownames(STAD_Integrate_celltype@meta.data) %in% rownames(DATA[which(DATA$NMF_Class %in% Com3),,drop=F]))]<-'Com3'
STAD_Integrate_celltype$NMF_Class[which(rownames(STAD_Integrate_celltype@meta.data) %in% rownames(DATA[which(DATA$NMF_Class %in% Com4),,drop=F]))]<-'Com4'

#DimPlot(STAD_Integrate_celltype,group.by="NMF_Class",reduction = "tsne",pt.size=0.1,label=TRUE)
Idents(STAD_Integrate_celltype)<-STAD_Integrate_celltype$NMF_Class

averger_exp <- AverageExpression(STAD_Integrate_celltype)[[1]]
averger_exp = averger_exp[,1:4]

genelist<-read.csv("/public/workspace/liuqzh/gastric_cancer/gene_list_all.csv",sep = ",",header = F)
genelist=genelist[which(!duplicated(genelist[,1])),]
rownames(genelist)=genelist[,1]

averger=averger_exp[which(rownames(averger_exp) %in% rownames(genelist)),]
sig_matrix=averger

#write.table(sig_matrix, "/public/workspace/liuqzh/gastric_cancer/cibersort/sig_matrix_Com_IIII.txt", row.names = T, sep = "\t")


rowscale <- results[,1:ncol(sig_matrix)]#只是相当于备份了一下results
rowscale <- rowscale[,apply(rowscale, 2, function(x){sum(x)>0})]#删除全是0的列
pheatmap(rowscale,
         scale = 'row',#按行标准化，不标准化就会按绝对值显示，很诡异
         cluster_col=T,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
         cluster_row=T,#是否对行聚类
         angle_col = "315")#调整X轴坐标的倾斜角度