bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ComplexHeatmap)

TCGA_STAD_Stemness_Score=read.table('/public/workspace/liuqzh/gastric_cancer/TCGA_subtype/TCGA_STAD_Stemness_Score.txt',header=T,sep='\t',row.names=1)
TCGA_STAD_Stemness_Score=TCGA_STAD_Stemness_Score[which(TCGA_STAD_Stemness_Score$sample_type == '1'),]
rownames(TCGA_STAD_Stemness_Score)=substr(rownames(TCGA_STAD_Stemness_Score),1,16)
TCGA_STAD_Stemness_Score$ID=rownames(TCGA_STAD_Stemness_Score)

MP1_MP2_MP3_TCGA=readRDS('/public/workspace/liuqzh/gastric_cancer/TCGA_subtype/MP1_MP2_MP3_TCGA_scrr.rds')
MP1_MP2_MP3_TCGA$ID=rownames(MP1_MP2_MP3_TCGA)
TCGA_STAD_Stemness_result=merge(TCGA_STAD_Stemness_Score,MP1_MP2_MP3_TCGA,'ID')

rownames(TCGA_STAD_Stemness_result)=TCGA_STAD_Stemness_result[,1]

data_f=TCGA_STAD_Stemness_result
p1=ggplot(data=data_f, aes(x=mRNAsi, y=MP_1))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x = mRNAsi, y = MP_1)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 0.5, label.y = 1) +  # 添加皮尔逊相关系数
	theme_minimal() +  # 使用简洁主题
	labs(
		title = "Relationship between MP1 score and mRNAsi",
		x = "mRNAsi",
		y = "MP_1"
	) +  # 添加图形标题和轴标题
	theme(
		plot.title = element_text(hjust = 0.5),  # 标题居中
		axis.text = element_text(color = "gray20"),  # 轴文字颜色
		axis.title = element_text(color = "gray20"),  # 轴标题颜色
		panel.grid.major = element_blank(),  # 移除主要网格线
		panel.grid.minor = element_blank(),  # 移除次要网格线
		panel.background = element_rect(fill = "white", color = "gray50")  # 背景颜色
	)

pdf('/public/workspace/liuqzh/gastric_cancer/TCGA_subtype/mRNAsi_MP1.pdf',height=7,width=7)
p1
dev.off()


