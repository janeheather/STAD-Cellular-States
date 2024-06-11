
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
setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/')
sur<-read.csv('GSE13861_YUHS_sur.csv',header=T,row.names=2)
gse<-read.table('GSE13861_matrix.txt',header=T,sep='\t',row.names=1)
ref<-read.table('GPL6884.txt',header=T,sep='\t',row.names=1)

GSE13861=gse[intersect(rownames(ref),rownames(gse)),]

ref$ID_REF=rownames(ref)
GSE13861$ID_REF=rownames(GSE13861)
GSE13861=merge(GSE13861,ref,'ID_REF')

final<-aggregate(GSE13861[,2:(ncol(GSE13861)-1)],list(GSE13861$Symbol),mean)
rownames(final)<-final[,1]
final<-final[,-1]
#final<-log2(final+1)
final=as.matrix(final)
sur<-arrange(sur,Patients_ID)
GSE13861_data=final

#saveRDS(GSE13861_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_data.rds')
GSE13861_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE13861_data<-GSE13861_data[rowSums(GSE13861_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE13861_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE13861_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE13861_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')
Result$GEO_ID=rownames(Result)
sur$GEO_ID=rownames(sur)

library(ggplot2)
library(ggpubr)
library(reshape2)
Result=merge(Result,sur,'GEO_ID')
Result$Recur=Result$Recurrence

Result$Recurrence='Na'
Result$Recurrence[Result$Recur==0]='Non_Recur'
Result$Recurrence[Result$Recur==1]='Recur'

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式

df_long <- melt(Result, measure.vars = c("MP1", "MP2", "MP3"), variable.name = "group", value.name = "Value")
##MP1
df=df_long[which(df_long$group == 'MP1'),]
df$group <- factor(df$Recurrence, levels = c("Recur", "Non_Recur"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = Value)) +
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
	labs(title = "GSE13861", x = "Group", y = "MP1 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Recur", "Non_Recur")), method='wilcox.test',
                            label = "p.signif", label.y = c(6))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_MP1.pdf',height=6,width=4)
print(p)
dev.off()

##MP2
df=df_long[which(df_long$group == 'MP2'),]
df$group <- factor(df$Recurrence, levels = c("Recur", "Non_Recur"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = Value)) +
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
	labs(title = "GSE13861", x = "Group", y = "MP2 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Recur", "Non_Recur")), method='wilcox.test',
                            label = "p.signif", label.y = c(2.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_MP2.pdf',height=6,width=4)
print(p)
dev.off()

##MP3
df=df_long[which(df_long$group == 'MP3'),]
df$group <- factor(df$Recurrence, levels = c("Recur", "Non_Recur"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = Value)) +
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
	labs(title = "GSE13861", x = "Group", y = "MP3 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Recur", "Non_Recur")), method='t.test',
                            label = "p.signif", label.y = c(3.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_MP3.pdf',height=6,width=4)
print(p)
dev.off()

#################
#################
#################
# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式

df_long <- melt(Result, measure.vars = c("MP1", "MP2", "MP3"), variable.name = "group", value.name = "Value")
df_long=df_long[which(df_long$Adjuvant_chemo == 'No'),]
##MP1
df=df_long[which(df_long$group == 'MP1'),]
df$group <- factor(df$Recurrence, levels = c("Recur", "Non_Recur"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = Value)) +
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
	labs(title = "GSE13861", x = "Group", y = "MP1 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Recur", "Non_Recur")), method='wilcox.test',
                            label = "p.signif", label.y = c(6))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_MP1_no_Adjuvant_chemo.pdf',height=6,width=4)
print(p)
dev.off()

##MP2
df=df_long[which(df_long$group == 'MP2'),]
df$group <- factor(df$Recurrence, levels = c("Recur", "Non_Recur"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = Value)) +
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
	labs(title = "GSE13861", x = "Group", y = "MP2 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Recur", "Non_Recur")), method='wilcox.test',
                            label = "p.signif", label.y = c(2.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_MP2_no_Adjuvant_chemo.pdf',height=6,width=4)
print(p)
dev.off()

##MP3
df=df_long[which(df_long$group == 'MP3'),]
df$group <- factor(df$Recurrence, levels = c("Recur", "Non_Recur"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = Value)) +
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
	labs(title = "GSE13861", x = "Group", y = "MP3 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Recur", "Non_Recur")), method='wilcox.test',
                            label = "p.signif", label.y = c(3.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_MP3_no_Adjuvant_chemo.pdf',height=6,width=4)
print(p)
dev.off()

