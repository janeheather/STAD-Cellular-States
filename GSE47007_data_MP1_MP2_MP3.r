
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
GSE47007_data=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE47007/GSE47007_series_matrix.txt",sep = "\t",header = T)
ref=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE47007/ref.txt",sep = "\t",header = T)
colnames(ref)[1]='ID_REF'
GSE47007_data<-merge(ref,GSE47007_data,by='ID_REF')

final<-aggregate(GSE47007_data[,3:ncol(GSE47007_data)],list(GSE47007_data[,2]),mean)

GSE47007_data=final[-c(1:14),]
rownames(GSE47007_data)=GSE47007_data[,1]
GSE47007_data=GSE47007_data[,-1]

#saveRDS(GSE47007_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE47007/GSE47007_data.rds')

GSE47007_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE47007/GSE47007_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE47007_data<-GSE47007_data[rowSums(GSE47007_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE47007_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE47007_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE47007_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)
Result$Subtype=''
Result$Subtype[grep('^Diffuse',rownames(Result))]='Diffuse'
Result$Subtype[grep('^Intestinal',rownames(Result))]='Intestinal'

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result

df_long <- melt(df, id.vars = c("Subtype"), measure.vars = c("MP1", "MP2", "MP3"), variable.name = "group", value.name = "Value")

df=df_long
df$group <- factor(df$group, levels = c("MP1", "MP2", "MP3"))  # 设置分组的顺序
df$Subtype <- factor(df$Subtype, levels = c("Diffuse","Intestinal"))  # 设置分组的顺序

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
	labs(title = "GSE47007", x = "Subtype", y = "MP1 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Diffuse","Intestinal")), 
                            method='wilcox.test',label = "p.signif", label.y = c(7))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE47007/MP1.pdf',height=6,width=3)
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
	labs(title = "GSE47007", x = "Subtype", y = "MP2 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Diffuse","Intestinal")), 
                            method='wilcox.test',label = "p.signif", label.y = c(3.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE47007/MP2.pdf',height=6,width=3)
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
	labs(title = "GSE47007", x = "Subtype", y = "MP3 Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("Diffuse","Intestinal")), 
                            method='wilcox.test',label = "p.signif", label.y = c(3.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE47007/MP3.pdf',height=6,width=3)
print(p)
dev.off()


