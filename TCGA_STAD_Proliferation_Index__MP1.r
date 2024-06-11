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

Proliferation_Index=read.table('/public/workspace/liuqzh/gastric_cancer/TCGA_subtype/Proliferation_Index.txt',header=T,sep='\t',row.names=1)

setwd('/public/workspace/liuqzh/gastric/survival/')
dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
#write.table(dat_exp,file="/public/workspace/liuqzh/gastric_cancer/bulk_data/TCGA/tcga_dat_expr.txt",row.names=T,col.names=T,quote=F,sep="\t")
sur<-read.table('/public/workspace/yumiao/STAD/tcga/TCGA-STAD.survival.tsv',header=T,sep='\t',row.names=1)

colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]
sur=sur[grep("-01A",rownames(sur)),]

sur=sur[intersect(rownames(sur),colnames(dat_exp)),]
dat_exp=dat_exp[,intersect(rownames(sur),colnames(dat_exp))]
dat_exp<-dat_exp[rowSums(dat_exp)>0,]

sur$OS.time=sur$OS.time/30
sur$OS[which(sur$OS.time >= 50)]=0

TPM=dat_exp
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/TCGA_subtype/Proliferation_Index_scr.csv",sep = ",",row.names = 1,header = F)
geneset<-c()
for (i in 1:nrow(genelist)) {
	aa<-genelist[i,]
	aa<-aa[!is.na(aa)]
	aa<-aa[!duplicated(aa)]
	aa<-aa[aa%in%rownames(TPM)]
	aa<-list(as.character(as.matrix(aa)))
	geneset<-c(geneset,aa)
}
names(geneset)<-rownames(genelist)
colCenter = function(m, by = 'mean') {
  m = as.matrix(m)
  if (by == 'mean')  by = T
  else if (by == 'median') by = matrixStats::colMedians(m)
  else stopifnot(is.numeric(by) & length(by) == ncol(m))
  scale(m, center = by, scale = F)
}

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)
colnames(ssGSEA_Score)='Proliferation_Index'
ssGSEA_Score$ID=rownames(ssGSEA_Score)
TCGA_STAD_Stemness_Score=ssGSEA_Score

MP1_MP2_MP3_TCGA=readRDS('/public/workspace/liuqzh/gastric_cancer/TCGA_subtype/MP1_MP2_MP3_TCGA_scrr.rds')
MP1_MP2_MP3_TCGA$ID=rownames(MP1_MP2_MP3_TCGA)
TCGA_STAD_Stemness_result=merge(TCGA_STAD_Stemness_Score,MP1_MP2_MP3_TCGA,'ID')

rownames(TCGA_STAD_Stemness_result)=TCGA_STAD_Stemness_result[,1]

data_f=TCGA_STAD_Stemness_result
p1=ggplot(data=data_f, aes(x=mRNAsi, y=MP_1))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x = Proliferation_Index, y = MP_1)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = -0.25, label.y = 1) +  # 添加皮尔逊相关系数
	theme_minimal() +  # 使用简洁主题
	labs(
		title = "MP1 score and Proliferation_Index",
		x = "Proliferation_Index",
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

pdf('/public/workspace/liuqzh/gastric_cancer/TCGA_subtype/Proliferation_Index_MP1.pdf',height=8,width=8)
p1
dev.off()


