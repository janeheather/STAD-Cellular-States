bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

GSE141657<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE141657/GSE141657_series_matrix.txt',header=T,sep='\t')
eff_length <- read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE141657/GPL23572-33493.txt',header=T,sep='\t')
colnames(eff_length)[1]='ID_REF'
GSE141657=merge(GSE141657,eff_length,'ID_REF')
final<-aggregate(GSE141657[,2:(ncol(GSE141657)-1)],list(GSE141657$GENE_SYMBOL),mean)
final=final[-c(1:28),]
rownames(final)<-final[,1]
final<-final[,-1]
GSE141657_data=final
saveRDS(GSE141657_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE141657/GSE141657_data.rds')
#########

GSE141657_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE141657/GSE141657_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE141657_data<-GSE141657_data[rowSums(GSE141657_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE141657_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE141657_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE141657_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')

library(ggplot2)
library(ggpubr)
library(reshape2)

Result$subtype=''
Result$subtype[grep('*Antrum',rownames(Result))]='Antrum'
Result$subtype[grep('*Corpus',rownames(Result))]='Corpus'
Result$subtype[grep('*Fundus',rownames(Result))]='Fundus'

Result$Color=''
Result$Color[grep('*_undiff',rownames(Result))]='undiff'
Result$Color[grep('*_diff',rownames(Result))]='diff'
# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result

df_long <- melt(df, measure.vars = c("MP1", "MP2", "MP3"), variable.name = "group", value.name = "Value")

######MP1
df=df_long[which(df_long$group == 'MP1'),]

df$subtype <- factor(df$subtype, levels = c("Fundus", "Corpus", "Antrum"))  # 设置分组的顺序
df$Color <- factor(df$Color, levels = c("diff", "undiff"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = subtype, y = Value, color = Color)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("diff" = "red", "undiff" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE141657", x = "Group", y = "MP1 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Fundus", "Corpus"), c("Corpus", "Antrum"), c("Fundus", "Antrum")), method='wilcox.test',
                            label = "p.signif", label.y = c(16.5, 17, 17.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE141657/GSE141657_MP1.pdf',height=6,width=4)
print(p)
dev.off()

######MP2
df=df_long[which(df_long$group == 'MP2'),]

df$subtype <- factor(df$subtype, levels = c("Fundus", "Corpus", "Antrum"))  # 设置分组的顺序
df$Color <- factor(df$Color, levels = c("diff", "undiff"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = subtype, y = Value, color = Color)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("diff" = "red", "undiff" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE141657", x = "Group", y = "MP2 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Fundus", "Corpus"), c("Corpus", "Antrum"), c("Fundus", "Antrum")), 
							method='wilcox.test',label = "p.signif", label.y = c(6, 6.5, 7))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE141657/GSE141657_MP2.pdf',height=6,width=4)
print(p)
dev.off()

######MP3
df=df_long[which(df_long$group == 'MP3'),]

df$subtype <- factor(df$subtype, levels = c("Fundus", "Corpus", "Antrum"))  # 设置分组的顺序
df$Color <- factor(df$Color, levels = c("diff", "undiff"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = subtype, y = Value, color = Color)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("diff" = "red", "undiff" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE141657", x = "Group", y = "MP3 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Fundus", "Corpus"), c("Corpus", "Antrum"), c("Fundus", "Antrum")), 
                            method='wilcox.test',label = "p.signif", label.y = c(9, 9.5, 10))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE141657/GSE141657_MP3.pdf',height=6,width=4)
print(p)
dev.off()
