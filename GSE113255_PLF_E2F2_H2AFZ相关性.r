
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
names_file=list.files('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE113255/')
file_path='/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE113255/'
file_all=NULL
for(i in names_file){
	print(i)
	tmp_file=read.table(paste(file_path,i,sep=''),header=T,sep='\t',row.names=1)
	tmp_file=as.matrix(tmp_file)
	file_all=cbind(file_all,tmp_file)
}
GSE113255_data=file_all[,-c(grep('N$',colnames(file_all)))]
colnames(GSE113255_data)=gsub('\\.','',colnames(GSE113255_data))

#saveRDS(GSE113255_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE113255/GSE113255_data.rds')
GSE113255_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE113255/GSE113255_data.rds')

dat_exp_H_E=GSE113255_data[which(rownames(GSE113255_data) %in% c('TUBA1B','E2F2')),]

dat_exp_H_E=as.data.frame(t(dat_exp_H_E))
data_f=dat_exp_H_E
p1=ggplot(data=data_f, aes(x=TUBA1B, y=E2F2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=TUBA1B, y=E2F2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 7.5, label.y = 6) +  # 添加皮尔逊相关系数
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

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE113255_TUBA1B_E2F2.pdf',height=5,width=5)
p1
dev.off()

##################

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
names_file=list.files('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE113255/')
file_path='/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE113255/'
file_all=NULL
for(i in names_file){
	print(i)
	tmp_file=read.table(paste(file_path,i,sep=''),header=T,sep='\t',row.names=1)
	tmp_file=as.matrix(tmp_file)
	file_all=cbind(file_all,tmp_file)
}
GSE113255_data=file_all[,-c(grep('N$',colnames(file_all)))]
colnames(GSE113255_data)=gsub('\\.','',colnames(GSE113255_data))

dat_exp_H_E=GSE113255_data[which(rownames(GSE113255_data) %in% c('H2AFZ','E2F2')),]

dat_exp_H_E=as.data.frame(t(dat_exp_H_E))
data_f=dat_exp_H_E
p1=ggplot(data=data_f, aes(x=H2AFZ, y=E2F2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=H2AFZ, y=E2F2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 6, label.y = 6) +  # 添加皮尔逊相关系数
	theme_minimal() +  # 使用简洁主题
	labs(
		title = "H2AFZ and E2F2 correlation",
		x = "H2AFZ",
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

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE113255_H2AFZ_E2F2.pdf',height=5,width=5)
p1
dev.off()




