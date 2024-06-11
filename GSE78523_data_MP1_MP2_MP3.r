
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(dplyr)
library(ggplot2)
library(ggtern)
library(GSVA)

set.seed(520)
GSE78523_data=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE78523/GSE78523_series_matrix.txt",sep = "\t",header = T)

ref=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE78523/GPL18990-3879.txt",sep = "\t",header = T)
colnames(ref)[1]='ID_REF'
GSE78523_data<-merge(ref,GSE78523_data,by='ID_REF')

final<-aggregate(GSE78523_data[,4:ncol(GSE78523_data)],list(GSE78523_data[,3]),mean)
final=final[-c(1:3),]

GSE78523_data=final
rownames(GSE78523_data)=GSE78523_data[,1]
GSE78523_data=GSE78523_data[,-1]

#saveRDS(GSE78523_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE78523/GSE78523_data.rds')

GSE78523_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE78523/GSE78523_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE78523_data<-GSE78523_data[rowSums(GSE78523_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE78523_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE78523_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE78523_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')

library(ggplot2)
library(ggpubr)
library(reshape2)
Result$Subtype=c()
Result$Subtype[grep('^Healthy_C',rownames(Result))]='NC'
Result$Subtype[grep('^IIM_C',rownames(Result))]='IIM_Con'
Result$Subtype[grep('^IIM_G',rownames(Result))]='IIM_GC'
Result$Subtype[grep('^CIM_C',rownames(Result))]='CIM_Con'
Result$Subtype[grep('^CIM_G',rownames(Result))]='CIM_GC'

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result

df_long <- melt(df, id.vars = c("Subtype"), measure.vars = c("MP1", "MP2", "MP3"), variable.name = "group", value.name = "Value")

df=df_long
df$group <- factor(df$group, levels = c("MP1", "MP2", "MP3"))  # 设置分组的顺序
df$Subtype <- factor(df$Subtype, levels = c("NC","IIM_Con","IIM_GC","CIM_Con","CIM_GC"))  # 设置分组的顺序

###MP1
MP1=df[which(df$group == 'MP1'),]
# 生成箱形图并添加散点图层
p <- ggplot(MP1, aes(x = Subtype, y = Value)) +
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
	labs(title = "GSE78523", x = "Subtype", y = "MP1 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("IIM_Con", "IIM_GC"), c("CIM_Con", "CIM_GC"),c("IIM_GC", "CIM_GC")), 
                            method='wilcox.test',label = "p.signif", label.y = c(2, 2,2.2))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE78523/MP1.pdf',height=6,width=8)
print(p)
dev.off()

###MP2
MP2=df[which(df$group == 'MP2'),]
# 生成箱形图并添加散点图层
p <- ggplot(MP2, aes(x = Subtype, y = Value)) +
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
	labs(title = "GSE78523", x = "Subtype", y = "MP2 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("IIM_Con", "IIM_GC"), c("CIM_Con", "CIM_GC"),c("IIM_GC", "CIM_GC")), 
                            method='wilcox.test',label = "p.signif", label.y = c(2.3, 2.3,2.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE78523/MP2.pdf',height=6,width=8)
print(p)
dev.off()

###MP3
MP3=df[which(df$group == 'MP3'),]
# 生成箱形图并添加散点图层
p <- ggplot(MP3, aes(x = Subtype, y = Value)) +
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
	labs(title = "GSE78523", x = "Subtype", y = "MP3 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("IIM_Con", "IIM_GC"), c("CIM_Con", "CIM_GC"),c("IIM_GC", "CIM_GC")), 
                            method='wilcox.test',label = "p.signif", label.y = c(2.5, 2.5,2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE78523/MP3.pdf',height=6,width=8)
print(p)
dev.off()


