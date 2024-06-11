
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

set.seed(520)
STAD_Tr_info<-readRDS("/public/workspace/liuqzh/gastric/STAD-PRJEB25780_info.Rds")
#STAD_Tr_info=as.data.frame(STAD_Tr_info[1:45,])
#rownames(STAD_Tr_info)=STAD_Tr_info[,1]

STAD_Tr_data<-readRDS("/public/workspace/liuqzh/gastric/STAD-PRJEB25780.Response.Rds")
STAD_Tr_data=as.data.frame(STAD_Tr_data)
rownames(STAD_Tr_data)=STAD_Tr_data[,1]
STAD_Tr_data=STAD_Tr_data[,-1]

genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
STAD_Tr_data<-STAD_Tr_data[rowSums(STAD_Tr_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(STAD_Tr_data),GENE_list, method='ssgsea')
GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(STAD_Tr_data),GENE_list, method='ssgsea')
GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(STAD_Tr_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')
Result$sample_id=rownames(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)
STAD_Tr_info=as.data.frame(STAD_Tr_info)
rownames(STAD_Tr_info)=STAD_Tr_info[,1]
STAD_Tr_info_Result=merge(STAD_Tr_info,Result,'sample_id')
rownames(STAD_Tr_info_Result)=STAD_Tr_info_Result[,1]

Result=STAD_Tr_info_Result
Result$Tumor_Normal='Tumor'
Result$Tumor_Normal[grep('-N$',Result$patient_name)]='Normal'

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result

df$Tumor_Normal <- factor(df$Tumor_Normal, levels = c("Tumor", "Normal"))  #设置分组的顺序
df$Sample=gsub('-N','',Result$patient_name)
#####
# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Tumor_Normal, y = MP1, color = Tumor_Normal)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") + #不显示离群值，且箱形图无填充
	scale_fill_manual() +  # 设置颜色
	geom_point(size=2,color='lightblue') +
	geom_line(aes(group = Sample), color = 'gray', lwd = 0.5) + #配对样本间连线
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "PRJEB25780_data", x = "Group", y = "MP1 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Tumor", "Normal")), method='wilcox.test',label = "p.signif")
pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/PRJEB25780/PRJEB25780_data_Tumor_Normal_MP1.pdf',height=6,width=3)
print(p)
dev.off()

#####
# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Tumor_Normal, y = MP2, color = Tumor_Normal)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	scale_fill_manual() +  # 设置颜色
	geom_point(size=2,color='lightblue') +
	geom_line(aes(group = Sample), color = 'gray', lwd = 0.5) +  #配对样本间连线
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "PRJEB25780_data", x = "Group", y = "MP2 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Tumor", "Normal")), method='wilcox.test',label = "p.signif")
pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/PRJEB25780/PRJEB25780_data_Tumor_Normal_MP2.pdf',height=6,width=3)
print(p)
dev.off()

#####
# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Tumor_Normal, y = MP3, color = Tumor_Normal)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	scale_fill_manual() +  # 设置颜色
	geom_point(size=2,color='lightblue') +
	geom_line(aes(group = Sample), color = 'gray', lwd = 0.5) +  #配对样本间连线
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "PRJEB25780_data", x = "Group", y = "MP3 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Tumor", "Normal")), method='wilcox.test',label = "p.signif")
pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/PRJEB25780/PRJEB25780_data_Tumor_Normal_MP3.pdf',height=6,width=3)
print(p)
dev.off()
###########
###########


