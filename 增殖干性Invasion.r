#############干性得分##############
{
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
colors=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")
dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(dat) <- "RNA"

data_single_cell_tumor=dat

gene <- c("DNMT3B","PFAS","XRCC5","HAUS6","TET1","IGF2BP1","PLAA",
"TEX10","MSH6","DLGAP5","SKIV2L2","SOHLH2","RRAS2","PAICS","CPSF3",
"LIN28B","IPO5","BMPR1A","ZNF788","ASCC3","FANCB","HMGA2","TRIM24",
"ORC1","HDAC2","HESX1","INHBE","MIS18A","DCUN1D5","MRPL3","CENPH",
"MYCN","HAUS1","GDF3","TBCE","RIOK2","BCKDHB","RAD1","NREP","ADH5",
"PLRG1","ROR1","RAB3B","DIAPH3","GNL2","FGF2","NMNAT2","KIF20A",
"CENPI","DDX1","XXYLT1","GPR176","BBS9","C14orf166","BOD1","CDC123",
"SNRPD3","FAM118B","DPH3","EIF2B3","RPF2","APLP1","DACT1","PDHB",
"C14orf119","DTD1","SAMM50","CCL26","MED20","UTP6","RARS2","ARMCX2",
"RARS","MTHFD2","DHX15","HTR7","MTHFD1L","ARMC9","XPOT","IARS","HDX",
"ACTRT3","ERCC2","TBC1D16","GARS","KIF7","UBE2K","SLC25A3","ICMT",
"UGGT2","ATP11C","SLC24A1","EIF2AK4","GPX8","ALX1","OSTC","TRPC4",
"HAS2","FZD2","TRNT1","MMADHC","SNX8","CDH6","HAT1","SEC11A","DIMT1","TM2D2","FST","GBE1")

DefaultAssay(data_single_cell_tumor) <- "RNA"
expr<- as.matrix(data_single_cell_tumor@assays$RNA@data)
expr<-expr[rowSums(expr)>0,]
genesets<-list(gene)
gsva.res<-gsva(expr,genesets,method = "ssgsea",parallel.sz = 10)
data_single_cell_tumor$stemnness_Score<-as.numeric(gsva.res)

test<-data_single_cell_tumor@meta.data
test$cluster<-factor(test$Subtype)

test=test[which(test$Subtype=="MP_1"|test$Subtype=="MP_2"|test$Subtype=="MP_3"),]
test$cluster<-factor(test$Subtype,levels=c("MP_1","MP_2","MP_3"))
Class <- list(c("MP_1","MP_2"),c("MP_2","MP_3"),c("MP_1","MP_3"))

colors=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")

# 生成箱形图并添加散点图层
p <- ggplot(test, aes(x=cluster,y=stemnness_Score,color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "Stemnness Score")  # 添加标签

pdf('/public/workspace/liuqzh/gastric_cancer/ste_score/stemnness_Score.pdf',width=6,height=5)
print(p)
dev.off()
}

#############增殖得分##############
{
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

dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(dat) <- "RNA"

data_single_cell_tumor=dat

gene <- c("AURKA","BUB1","CCNB1","CCND1","CCNE1","DEK",
		  "E2F1","FEN1","FOXM1","H2AFZ","HMGB2","MCM2",
		  "MCM3","MCM4","MCM5","MCM6","MKI67","MYBL2",
		  "PCNA","PLK1","TOP2A","TYMS","ZWINT")

DefaultAssay(data_single_cell_tumor) <- "RNA"
expr<- as.matrix(data_single_cell_tumor@assays$RNA@data)
expr<-expr[rowSums(expr)>0,]
genesets<-list(gene)
gsva.res<-gsva(expr,genesets,method = "ssgsea",parallel.sz = 10)
data_single_cell_tumor$stemnness_Score<-as.numeric(gsva.res)

test<-data_single_cell_tumor@meta.data         
test$cluster<-factor(test$Subtype)

test=test[which(test$Subtype=="MP_1"|test$Subtype=="MP_2"|test$Subtype=="MP_3"),]
test$cluster<-factor(test$Subtype,levels=c("MP_1","MP_2","MP_3"))
Class <- list(c("MP_1","MP_2"),c("MP_2","MP_3"),c("MP_1","MP_3"))

colors=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")

# 生成箱形图并添加散点图层 
p <- ggplot(test, aes(x=cluster,y=stemnness_Score,color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "Proliferation Score")  # 添加标签

pdf('/public/workspace/liuqzh/gastric_cancer/ste_score/Proliferation_Score.pdf',width=6,height=5)
print(p)
dev.off()
}

#############Invasion##############
{
bytlib load R-4.0.2
R

library(Seurat)
library(monocle)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)

dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(dat) <- "RNA"

data_single_cell_tumor=dat

gene <- c('TGFB1','CSF1','EGF','VEGFA','PTGS2','EREG','ANGPT2','MMP1','MMP2','MMP3','MMP10')

DefaultAssay(data_single_cell_tumor) <- "RNA"
expr<- as.matrix(data_single_cell_tumor@assays$RNA@data)
expr<-expr[rowSums(expr)>0,]
genesets<-list(gene)
gsva.res<-gsva(expr,genesets,method = "ssgsea",parallel.sz = 10)
data_single_cell_tumor$stemnness_Score<-as.numeric(gsva.res)

test<-data_single_cell_tumor@meta.data         
test$cluster<-factor(test$Subtype)

test=test[which(test$Subtype=="MP_1"|test$Subtype=="MP_2"|test$Subtype=="MP_3"),]
test$cluster<-factor(test$Subtype,levels=c("MP_1","MP_2","MP_3"))
Class <- list(c("MP_1","MP_2"),c("MP_2","MP_3"),c("MP_1","MP_3"))

colors=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")

# 生成箱形图并添加散点图层
p <- ggplot(test, aes(x=cluster,y=stemnness_Score,color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "Invasion Score")  # 添加标签

pdf('/public/workspace/liuqzh/gastric_cancer/ste_score/Invasion_score.pdf',width=6,height=5)
print(p)
dev.off()
}



