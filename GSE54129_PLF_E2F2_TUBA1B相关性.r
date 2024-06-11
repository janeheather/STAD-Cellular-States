
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(ggplot2)
library(ggtern)
library(GSVA)
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(estimate)

set.seed(520)
GSE54129_data=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE54129/GSE54129_matrix.txt",sep = "\t",header = T)
GPL570_ref=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE54129/GPL570.txt",sep = "\t",header = T)
colnames(GPL570_ref)[1]='ID_REF'

GSE54129_data=merge(GPL570_ref,GSE54129_data,'ID_REF')
GSE54129_data=GSE54129_data[,-1]

final<-aggregate(GSE54129_data[,2:ncol(GSE54129_data)],list(GSE54129_data[,1]),mean)
final=final[-1,]
rownames(final)=final[,1]
final=final[,-1]
GSE54129_data=final

#saveRDS(GSE54129_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE54129/GSE54129_data.rds')
GSE54129_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE54129/GSE54129_data.rds')

GSE54129_data_normal=GSE54129_data[,grep('^normal_',colnames(GSE54129_data))]
GSE54129_data_tumor=GSE54129_data[,grep('^tumor_',colnames(GSE54129_data))]

dat_exp_H_E=GSE54129_data_tumor[which(rownames(GSE54129_data_tumor) %in% c('TUBA1B','E2F2')),]

dat_exp_H_E=as.data.frame(t(dat_exp_H_E))
data_f=dat_exp_H_E
p1=ggplot(data=data_f, aes(x=TUBA1B, y=E2F2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=TUBA1B, y=E2F2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 10.5, label.y = 5.5) +  # 添加皮尔逊相关系数
	theme_minimal() +  # 使用简洁主题
	labs(
		title = "TUBA1B and E2F2 correlation",
		x = "TUBA1B",
		y = "E2F2"
	) +  # 添加图形标题和轴标题
	theme(
		plot.title = element_text(hjust = 0.5),  # 标题居中
		axis.text = element_text(color = "gray20"),  # 轴文字颜色
		axis.title = element_text(color = "gray20"),  # 轴标题颜色
		panel.grid.major = element_blank(),  # 移除主要网格线
		panel.grid.minor = element_blank(),  # 移除次要网格线
		panel.background = element_rect(fill = "white", color = "gray50")  # 背景颜色
	)
p1

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE54129_TUBA1B_E2F2.pdf',height=5,width=5)
p1
dev.off()
