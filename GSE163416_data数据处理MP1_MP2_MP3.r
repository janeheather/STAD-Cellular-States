bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

GSE163416_CAG_IM<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE163416/GSE163416_N-CAG+IM.all.mRNA.txt',header=T,sep='\t')
GSE163416_Dys<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE163416/GSE163416_N-Dys.all.mRNA.txt',header=T,sep='\t')
GSE163416_GC<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE163416/GSE163416_N-GC.all.mRNA.txt',header=T,sep='\t')

GSE163416_CAG_IM<-aggregate(GSE163416_CAG_IM[,2:ncol(GSE163416_CAG_IM)],list(GSE163416_CAG_IM[,1]),mean)
GSE163416_Dys<-aggregate(GSE163416_Dys[,2:ncol(GSE163416_Dys)],list(GSE163416_Dys[,1]),mean)
GSE163416_GC<-aggregate(GSE163416_GC[,2:ncol(GSE163416_GC)],list(GSE163416_GC[,1]),mean)

GSE163416_data=merge(GSE163416_CAG_IM,GSE163416_Dys,'Group.1')
GSE163416_data=merge(GSE163416_data,GSE163416_GC,'Group.1')

GSE163416_data=GSE163416_data[,-c(2,3,4,8,9,10,14,15,16)]
rownames(GSE163416_data)=GSE163416_data[,1]
GSE163416_data=GSE163416_data[,-1]

colnames(GSE163416_data)=gsub('G.I','G_I',colnames(GSE163416_data))
GSE163416_data=GSE163416_data[-c(1,2),]
saveRDS(GSE163416_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE163416/GSE163416_data.rds')

set.seed(520)
GSE163416_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE163416/GSE163416_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE163416_data<-GSE163416_data[rowSums(GSE163416_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE163416_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE163416_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE163416_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')

library(ggplot2)
library(ggpubr)
library(reshape2)

Result$subtype=c()
Result$subtype[grep('^CAG',rownames(Result))]='CAG_IM'
Result$subtype[grep('^Dys',rownames(Result))]='Dys'
Result$subtype[grep('^GC',rownames(Result))]='GC'

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result
df_long <- melt(df, measure.vars = c("MP1", "MP2", "MP3"), variable.name = "group", value.name = "Value")

###MP1
df=df_long[which(df_long$group == 'MP1'),]
df$subtype <- factor(df$subtype, levels = c("CAG_IM", "Dys", "GC"))  # 设置分组的顺序
# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = subtype, y = Value)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("red" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE163416", x = "Group", y = "MP1 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("CAG_IM", "Dys"), c("Dys", "GC"), c("CAG_IM", "GC")), method='wilcox.test',
                            label = "p.signif", label.y = c(3.5, 4, 4.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE163416/GSE163416_MP1.pdf',height=6,width=4)
print(p)
dev.off()

#######
###MP2
df=df_long[which(df_long$group == 'MP2'),]
df$subtype <- factor(df$subtype, levels = c("CAG_IM", "Dys", "GC"))  # 设置分组的顺序
# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = subtype, y = Value)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("red" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE163416", x = "Group", y = "MP2 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("CAG_IM", "Dys"), c("Dys", "GC"), c("CAG_IM", "GC")), method='wilcox.test',
                            label = "p.signif", label.y = c(3, 3.5, 4))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE163416/GSE163416_MP2.pdf',height=6,width=4)
print(p)
dev.off()

###MP3
df=df_long[which(df_long$group == 'MP3'),]
df$subtype <- factor(df$subtype, levels = c("CAG_IM", "Dys", "GC"))  # 设置分组的顺序
# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = subtype, y = Value)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("red" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE163416", x = "Group", y = "MP3 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("CAG_IM", "Dys"), c("Dys", "GC"), c("CAG_IM", "GC")), method='wilcox.test',
                            label = "p.signif", label.y = c(4.5, 5, 5.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE163416/GSE163416_MP3.pdf',height=6,width=4)
print(p)
dev.off()
