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

setwd("/public/workspace/liuqzh/gastric_cancer/bulk_data/TCGA/")
##将准备好的表达谱保存为txt格式，这里是用ncbiid，如果是用genesymbol,改成id="GeneSymbol"即可
#filterCommonGenes(input.f="tcga_dat_expr.txt", output.f="tcga.gct", id="GeneSymbol")
#estimateScore(input.ds="tcga.gct", output.ds="tcga_estimate_score.gct", platform="affymetrix")
#estimate_score <- read.table("tcga_estimate_score.gct", skip = 2, header = TRUE)
#write.csv(estimate_score,"tcga_est.csv",row.names = FALSE)

set.seed(520)
TPM=dat_exp

genelist=c('PDGFRB','BGN')
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)

data_f=ssGSEA_Score
p1=ggplot(data=data_f, aes(x=PDGFRB, y=BGN))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x = PDGFRB, y = BGN)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 2, label.y = 10) +  # 添加皮尔逊相关系数
	theme_minimal() +  # 使用简洁主题
	labs(
		title = "BGN and PDGFRB correlation",
		x = "PDGFRB",
		y = "BGN"
	) +  # 添加图形标题和轴标题
	theme(
		plot.title = element_text(hjust = 0.5),  # 标题居中
		axis.text = element_text(color = "gray20"),  # 轴文字颜色
		axis.title = element_text(color = "gray20"),  # 轴标题颜色
		panel.grid.major = element_blank(),  # 移除主要网格线
		panel.grid.minor = element_blank(),  # 移除次要网格线
		panel.background = element_rect(fill = "white", color = "gray50")  # 背景颜色
	)

pdf('/public/workspace/liuqzh/gastric_cancer/figure3/BGN_and_PDGFRB_correlation.pdf',height=7,width=7)
p1
dev.off()


