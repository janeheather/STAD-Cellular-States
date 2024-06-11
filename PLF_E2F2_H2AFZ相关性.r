####correlation####
###########################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)

setwd('/public/workspace/liuqzh/gastric/survival/')
dat_exp<-read.table('/public/workspace/liuqzh/gastric_cancer/cibersort/STAD_FPKM.txt',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

dat_exp_H_E=dat_exp[which(rownames(dat_exp) %in% c('H2AFZ','E2F2')),]
dat_exp_H_E=as.data.frame(t(dat_exp_H_E))
data_f=dat_exp_H_E
p1=ggplot(data=data_f, aes(x=H2AFZ, y=E2F2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=H2AFZ, y=E2F2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 4, label.y = 4) +  # 添加皮尔逊相关系数
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

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_TCGA_STAD_H2AFZ_E2F2.pdf',height=5,width=5)
p1
dev.off()
}

###########################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

setwd('/public/workspace/yumiao/STAD/bulk/GSE62254/')
sur<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE62254/GSE62254.csv',header=T)
gse<-read.table('GSE62254_series_matrix (1).txt',header=T,sep='\t')
ref<-read.table('GPL570-55999 (1).txt',header=T,sep='\t')

colnames(ref)[1]<-'ID_REF'
GSE62254<-merge(ref,gse,by='ID_REF')
final<-aggregate(GSE62254[,3:ncol(GSE62254)],list(GSE62254[,2]),mean)
final<-final[-1,]
rownames(final)<-final[,1]
final<-final[,-1]
#final<-log2(final+1)
sur<-arrange(sur,GEO_ID)
dat_exp=final

set.seed(12345)

dat_exp_H_E=dat_exp[which(rownames(dat_exp) %in% c('H2AFZ','E2F2')),]
dat_exp_H_E=as.data.frame(t(dat_exp_H_E))
data_f=dat_exp_H_E
p1=ggplot(data=data_f, aes(x=H2AFZ, y=E2F2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=H2AFZ, y=E2F2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 3.2, label.y = 2) +  # 添加皮尔逊相关系数
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

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE62254_H2AFZ_E2F2.pdf',height=5,width=5)
p1
dev.off()
}

###########################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE26253/')
sur<-read.csv('GSE26253_SMC_sur.csv',header=T,row.names=2)
gse<-read.table('GSE26253_matrix.txt',header=T,sep='\t')
ref<-read.table('GPL8432-11703.txt',header=T,sep='\t')

colnames(ref)[1]<-'ID_REF'
GSE26253<-merge(ref,gse,by='ID_REF')
final<-aggregate(GSE26253[,3:ncol(GSE26253)],list(GSE26253$Symbol),mean)
rownames(final)<-final[,1]
final<-final[,-1]
#final<-log2(final+1)
final=as.matrix(final)
sur<-arrange(sur,Patients_ID)

set.seed(12345)
dat_exp=final

dat_exp_H_E=dat_exp[which(rownames(dat_exp) %in% c('H2AFZ','E2F2')),]
dat_exp_H_E=as.data.frame(t(dat_exp_H_E))
data_f=dat_exp_H_E
p1=ggplot(data=data_f, aes(x=H2AFZ, y=E2F2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=H2AFZ, y=E2F2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 3.2, label.y = 2) +  # 添加皮尔逊相关系数
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

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE26253_H2AFZ_E2F2.pdf',height=5,width=5)
p1
dev.off()
}

###########################
{
bytlib load gcc
bytlib load R-4.0.2
R

library(GSVA)
library(survival)
library(survminer)
library(dplyr)
 
cli<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE84426_infor.csv',header=T)
expr<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE84426_matrix.csv',header=T)
anno<-read.table('/public/workspace/yumiao/STAD/bulk/GPL6947-13512.txt',sep='\t',header=T)
colnames(anno)[1]<-'ID_REF'
dat<-merge(anno,expr,by='ID_REF')
final<-aggregate(dat[,3:ncol(dat)],list(dat$Symbol),mean)
final<-final[-1,]
rownames(final)<-final[,1]
final<-final[,-1]

final<-log2(final+1)

set.seed(12345)
dat_exp=final

dat_exp_H_E=dat_exp[which(rownames(dat_exp) %in% c('H2AFZ','E2F2')),]
dat_exp_H_E=as.data.frame(t(dat_exp_H_E))
data_f=dat_exp_H_E
p1=ggplot(data=data_f, aes(x=H2AFZ, y=E2F2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=H2AFZ, y=E2F2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 3.2, label.y = 2) +  # 添加皮尔逊相关系数
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

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE84426_H2AFZ_E2F2.pdf',height=5,width=5)
p1
dev.off()
}

#############
{
bytlib load R-4.0.2
bytlib load gcc
R

library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

GSE15459_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE15459/GSE15459_data.rds')
GSE15459_info=read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE15459/GSE15459_data_clinical_information.csv',row.names=1)
GSE15459_info$sample=rownames(GSE15459_info)
cli=GSE15459_info

set.seed(12345)
dat_exp=GSE15459_data
dat_exp<-dat_exp[rowSums(dat_exp)>0,]

dat_exp_H_E=dat_exp[which(rownames(dat_exp) %in% c('H2AFZ','E2F2')),]
dat_exp_H_E=as.data.frame(t(dat_exp_H_E))
data_f=dat_exp_H_E
p1=ggplot(data=data_f, aes(x=H2AFZ, y=E2F2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=H2AFZ, y=E2F2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 12, label.y = 9) +  # 添加皮尔逊相关系数
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

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE15459_H2AFZ_E2F2.pdf',height=5,width=5)
p1
dev.off()
}

