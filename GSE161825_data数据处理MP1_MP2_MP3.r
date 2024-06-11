
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
GSE161825_data=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE161825/GSE161825_data.txt",sep = "\t",header = T)
final<-aggregate(GSE161825_data[,2:ncol(GSE161825_data)],list(GSE161825_data[,1]),mean)
GSE161825_data=final
rownames(GSE161825_data)=GSE161825_data[,1]
GSE161825_data=GSE161825_data[,-1]
GSE161825_data=log2(GSE161825_data+1)
GSE161825_data=GSE161825_data[-c(1:18,19:23,24:27),]

#saveRDS(GSE161825_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE161825/GSE161825_data.rds')
GSE161825_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE161825/GSE161825_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE161825_data<-GSE161825_data[rowSums(GSE161825_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE161825_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE161825_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE161825_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')

library(ggplot2)
library(ggpubr)
library(reshape2)
Result$Subtype=''
Result$Subtype[grep('^GCSC.Con',rownames(Result))]='Control'
Result$Subtype[grep('^GCSC.Tip.',rownames(Result))]='Tip'

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result

df_long <- melt(df, id.vars = c("Subtype"), measure.vars = c("MP1", "MP2", "MP3"), variable.name = "group", value.name = "Value")
df_long$Value=scale(df_long$Value)

df=df_long
df$group <- factor(df$group, levels = c("MP1", "MP2", "MP3"))  # 设置分组的顺序
df$Subtype <- factor(df$Subtype, levels = c("Control","Tip"))  # 设置分组的顺序

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
	labs(title = "GSE161825", x = "Group", y = "MP Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP1", "MP2"), c("MP2", "MP3"), c("MP1", "MP3")), 
                            label = "p.signif", label.y = c(1.5, 2, 2.5))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE161825/GSE161825_MP1_MP2_MP3.pdf',height=6,width=4)
print(p)
dev.off()

