
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
#读取scRNA
dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(dat) <- "RNA"

Single_RNA=as.matrix(dat@assays$RNA@data)
Single_RNA<-Single_RNA[rowSums(Single_RNA)>0,]

Hallmarks<-readRDS('/public/workspace/liuqzh/gastric_cancer/geneList_msigdb_Hallmarks.RDS')

###HALLMARK_HYPOXIA
{
GENE_list=list(Hallmarks$HALLMARK_HYPOXIA)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_HYPOXIA=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))
 
# 生成箱形图并添加散点图层
p <-ggplot(df, aes(x = Subtype, y = HALLMARK_HYPOXIA, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK HYPOXIA Score")  # 添加标签

# 计算并添加指定组间的p值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_HYPOXIA.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
{
GENE_list=list(Hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK EMT Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_EMT.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_TNFA_SIGNALING_VIA_NFKB
{
GENE_list=list(Hallmarks$HALLMARK_TNFA_SIGNALING_VIA_NFKB)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_TNFA_SIGNALING_VIA_NFKB=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_TNFA_SIGNALING_VIA_NFKB, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK TNFA SIGNALING VIA NFKB Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_TNFA.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_INFLAMMATORY_RESPONSE
{
GENE_list=list(Hallmarks$HALLMARK_INFLAMMATORY_RESPONSE)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_INFLAMMATORY_RESPONSE=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_INFLAMMATORY_RESPONSE, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK INFLAMMATORY RESPONSE Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_INFLAMMATORY_RESPONSE.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_PI3K_AKT_MTOR_SIGNALING
{
GENE_list=list(Hallmarks$HALLMARK_PI3K_AKT_MTOR_SIGNALING)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_PI3K_AKT_MTOR_SIGNALING=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_PI3K_AKT_MTOR_SIGNALING, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK PI3K AKT MTOR SIGNALING Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_PI3K_AKT_MTOR_SIGNALING.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_WNT_BETA_CATENIN_SIGNALING
{
GENE_list=list(Hallmarks$HALLMARK_WNT_BETA_CATENIN_SIGNALING)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_WNT_BETA_CATENIN_SIGNALING=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_WNT_BETA_CATENIN_SIGNALING, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK WNT BETA CATENIN SIGNALING Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_WNT_BETA_CATENIN_SIGNALING.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_E2F_TARGETS
{
GENE_list=list(Hallmarks$HALLMARK_E2F_TARGETS)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_E2F_TARGETS=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_E2F_TARGETS, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK E2F TARGETS Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_E2F_TARGETS.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_P53_PATHWAY
{
GENE_list=list(Hallmarks$HALLMARK_P53_PATHWAY)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_P53_PATHWAY=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_P53_PATHWAY, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK P53 PATHWAY Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_P53_PATHWAY.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_APICAL_JUNCTION
{
GENE_list=list(Hallmarks$HALLMARK_APICAL_JUNCTION)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_APICAL_JUNCTION=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_APICAL_JUNCTION, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK APICAL JUNCTION Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_APICAL_JUNCTION.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_ANGIOGENESIS
{
GENE_list=list(Hallmarks$HALLMARK_ANGIOGENESIS)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_ANGIOGENESIS=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_ANGIOGENESIS, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK ANGIOGENESIS Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_ANGIOGENESIS.pdf',height=6,width=5)
print(p)
dev.off()
}

###HALLMARK_NOTCH_SIGNALING
{
GENE_list=list(Hallmarks$HALLMARK_NOTCH_SIGNALING)
score=t(gsva(as.matrix(Single_RNA),GENE_list, method='ssgsea'))

dat@meta.data$HALLMARK_NOTCH_SIGNALING=score

Result=dat@meta.data[,16:17]
Result=as.data.frame(Result)

library(ggplot2)
library(ggpubr)
library(reshape2)

df=Result
df$Subtype=factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = HALLMARK_NOTCH_SIGNALING, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	#geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("MP_1" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA_seq", x = "Group", y = "HALLMARK NOTCH SIGNALING Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"), c("MP_2", "MP_3"), c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 2.6, 2.7))

pdf('/public/workspace/liuqzh/gastric_cancer/ssgsea/HALLMARK_NOTCH_SIGNALING.pdf',height=6,width=5)
print(p)
dev.off()
}