
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
GSE198136_data=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE198136/GSE198136.txt",sep = "\t",header = T)
rownames(GSE198136_data)=GSE198136_data$gene

GSE198136_data=GSE198136_data[,-1]
GSE198136_data=log2(GSE198136_data+1)
#saveRDS(GSE198136_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE198136/GSE198136_data.rds')

GSE198136_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE198136/GSE198136_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE198136_data<-GSE198136_data[rowSums(GSE198136_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE198136_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE198136_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE198136_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')

library(ggplot2)
library(ggpubr)
library(reshape2)

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result

df$Treat[grep('d1$',rownames(df))]='D1'
df$Treat[grep('d15$',rownames(df))]='D15'

df$Treat <- factor(df$Treat, levels = c("D1", "D15"))  # 设置分组的顺序

df$Sample='Other'
df$Sample[grep('^X1023',rownames(df))]='X1023'
df$Sample[grep('^X1024',rownames(df))]='X1024'
df$Sample[grep('^X1027',rownames(df))]='X1027'
df$Sample[grep('^X1028',rownames(df))]='X1028'
df$Sample[grep('^X1029',rownames(df))]='X1029'
df$Sample[grep('^X1031',rownames(df))]='X1031'
df$Sample[grep('^X1032',rownames(df))]='X1032'
df$Sample[grep('^X1034',rownames(df))]='X1034'
df$Sample[grep('^X1035',rownames(df))]='X1035'
df$Sample[grep('^X1038',rownames(df))]='X1038'
df$Sample[grep('^X1040',rownames(df))]='X1040'
df$Sample[grep('^X1042',rownames(df))]='X1042'
df$Sample[grep('^X1043',rownames(df))]='X1043'
df$Sample[grep('^X1046',rownames(df))]='X1046'
df$Sample[grep('^X1047',rownames(df))]='X1047'
df$Sample[grep('^X1048',rownames(df))]='X1048'
df$Sample[grep('^X1049',rownames(df))]='X1049'
df$Sample[grep('^X1050',rownames(df))]='X1050'
df$Sample[grep('^X1052',rownames(df))]='X1052'
df$Sample[grep('^X1054',rownames(df))]='X1054'
df$Sample[grep('^X2005',rownames(df))]='X2005'
df$Sample[grep('^X2006',rownames(df))]='X2006'
df$Sample[grep('^X2007',rownames(df))]='X2007'
# 生成箱形图并添加散点图层
df=as.data.frame(df)

df=df[which(df$Sample %in% c('X1027','X1028','X1031','X1032',
							 'X1034','X1035','X1038','X1047',
							 'X1049','X1050','X1052')),]
#####MP1
p <- ggplot(df, aes(x = Treat, y = MP1, color = Treat)) +
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
	labs(title = "GSE198136_data", x = "Group", y = "MP1 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("D1", "D15")),
							method='wilcox.test',paired=TRUE,label = "p.signif")

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE198136/GSE198136_D1_D15_MP1.pdf',height=6,width=3)
print(p)
dev.off()

#####MP2
p <- ggplot(df, aes(x = Treat, y = MP2, color = Treat)) +
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
	labs(title = "GSE198136_data", x = "Group", y = "MP2 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("D1", "D15")),
							method='wilcox.test',paired=TRUE,label = "p.signif")

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE198136/GSE198136_D1_D15_MP2.pdf',height=6,width=3)
print(p)
dev.off()

#####MP3
p <- ggplot(df, aes(x = Treat, y = MP3, color = Treat)) +
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
	labs(title = "GSE198136_data", x = "Group", y = "MP3 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("D1", "D15")),
							method='wilcox.test',paired=TRUE,label = "p.signif")

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE198136/GSE198136_D1_D15_MP3.pdf',height=6,width=3)
print(p)
dev.off()