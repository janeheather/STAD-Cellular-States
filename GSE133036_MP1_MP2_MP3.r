
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
GSE133036_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE133036/GSE133036_data.rds')
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv",sep = ",",header = T)
GSE133036_data<-GSE133036_data[rowSums(GSE133036_data)>0,]

GENE_list=list(genelist$MP_1)
ssGSEA_Score_MP1=gsva(as.matrix(GSE133036_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_2)
ssGSEA_Score_MP2=gsva(as.matrix(GSE133036_data),GENE_list, method='ssgsea')

GENE_list=list(genelist$MP_3)
ssGSEA_Score_MP3=gsva(as.matrix(GSE133036_data),GENE_list, method='ssgsea')

Result=rbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
Result=as.data.frame(t(Result))
colnames(Result)=c('MP1','MP2','MP3')

library(ggplot2)
library(ggpubr)
library(reshape2)
Result$AQP=c()
Result$AQP[grep('P$',rownames(Result))]='AQP_P'
Result$AQP[grep('N$',rownames(Result))]='AQP_N'
# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result
df$Color <- ifelse(df$AQP == 'AQP_P', "red", "lightblue")  # 设置颜色条件

df_long <- melt(df, id.vars = c("AQP", "Color"), measure.vars = c("MP1", "MP2", "MP3"), variable.name = "group", value.name = "Value")
df=df_long
df$group <- factor(df$group, levels = c("MP1", "MP2", "MP3"))  # 设置分组的顺序
df$AQP <- factor(df$AQP, levels = c("AQP_N", "AQP_P"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = Value, color = Color)) +
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
	labs(title = "GSE133036", x = "Group", y = "MP Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP1", "MP2"), c("MP2", "MP3"), c("MP1", "MP3")), 
                            label = "p.signif", label.y = c(13, 13.5, 14))

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE133036/MP1_MP2_MP3.pdf',height=6,width=4)
print(p)
dev.off()
