
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
GSE117420_data=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE117420/GSE117420_data.txt",sep = "\t",header = T)
GPL10787_ref=read.table("/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE117420/GPL10787-9758.txt",sep = "\t",header = T)
colnames(GPL10787_ref)[1]='ID_REF'

GSE117420_data=merge(GPL10787_ref,GSE117420_data,'ID_REF')
GSE117420_data=GSE117420_data[,-1]

final<-aggregate(GSE117420_data[,2:ncol(GSE117420_data)],list(GSE117420_data[,1]),mean)
final=final[-1,]
rownames(final)=final[,1]
final=final[,-1]
GSE117420_data=final

#saveRDS(GSE117420_data,'/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE117420/GSE117420_data.rds')
GSE117420_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE117420/GSE117420_data.rds')

Result=GSE117420_data[which(rownames(GSE117420_data) %in% c('Tuba1b')),]

Result=as.data.frame(t(Result))

library(ggplot2)
library(ggpubr)
library(reshape2)

Result$Subtype='Other'
Result$Subtype[grep('^Liver_WT',rownames(Result))]='E2F2 WT'
Result$Subtype[grep('^Liver_E2F2',rownames(Result))]='E2F2 KO'

df=as.data.frame(Result)

p <- ggplot(df, aes(x = Subtype, y = Tuba1b, color = Subtype)) +
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
	labs(title = "GSE117420_data", x = "Group", y = "Expression of Tuba1b")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("E2F2 KO", "E2F2 WT")), method='wilcox.test',
                            label = "p.signif")

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE117420/GSE117420_E2F2_tuba1b.pdf',height=6,width=3)
print(p)
dev.off()
#################

