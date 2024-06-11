bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(reshape2)
library(dplyr)
library(pheatmap)
library(philentropy)
library(tibble)
library(tidyr)
library(patchwork)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(stringr)
library("Rtsne")
library(RColorBrewer)
library(scales)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ggplot2)
library(ggpubr)
library(rlang)
source('/public/workspace/liuqzh/ssgseaMOD.r')

setwd('/public/workspace/liuqzh/ssgsea_tmp/mp1/')
dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
colnames(dat_exp)=gsub("\\.","-",colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

marker<-read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv',header=T,row.names=1)
Gene=marker$Maligant_MP_3

exp=dat_exp
mod.generate(Gene,out=paste0('/public/workspace/liuqzh/ssgsea_tmp/mp3/','Maligant_MP_3','.mod'))
mod <- mod.analyze2(exp,'Maligant_MP_3','/public/workspace/liuqzh/ssgsea_tmp/mp3/',permN=100000)
rownames(mod)=gsub("\\.","-",rownames(mod))

#write.table(mod,file="/public/workspace/liuqzh/gastric_cancer/Maligant_MP_3_TCGA.txt",row.names=T,col.names=T,quote=F,sep="\t")
##
marker<-read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv',header=T,row.names=1)
Gene=marker$Maligant_MP_2

exp=dat_exp
mod.generate(Gene,out=paste0('/public/workspace/liuqzh/ssgsea_tmp/mp2/','Maligant_MP_2','.mod'))
mod <- mod.analyze2(exp,'Maligant_MP_2','/public/workspace/liuqzh/ssgsea_tmp/mp2/',permN=100000)
rownames(mod)=gsub("\\.","-",rownames(mod))
##
marker<-read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv',header=T,row.names=1)
Gene=marker$Maligant_MP_1

exp=dat_exp
mod.generate(Gene,out=paste0('/public/workspace/liuqzh/ssgsea_tmp/mp1/','Maligant_MP_1','.mod'))
mod <- mod.analyze2(exp,'Maligant_MP_1','/public/workspace/liuqzh/ssgsea_tmp/mp1/',permN=100000)
rownames(mod)=gsub("\\.","-",rownames(mod))

##################
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(reshape2)
library(dplyr)
library(pheatmap)
library(philentropy)
library(tibble)
library(tidyr)
library(patchwork)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(stringr)
library("Rtsne")
library(RColorBrewer)
library(scales)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ggplot2)
library(ggpubr)
library(rlang)

TCGA=read.table('/public/workspace/liuqzh/gastric_cancer/TCGA_MP_subtype_SCORE.txt',header=T,sep='\t',row.names=1)

ha = HeatmapAnnotation(Subtype=as.character(TCGA$Subtype), col = list(Subtype= c("MP2" = '#509285', "MP1" = '#BA3630',"MP3" = '#A2855C')))

Heatmap(t(TCGA[,c(25:27)]),cluster_columns = F,cluster_rows = F,show_column_names = F,name = "heat",top_annotation = ha,col = c('#5A8FCA','#F2F2F0','#E31A1C'))

data_TCGA=read.table('/public/workspace/liuqzh/gastric_cancer/info_data/TCGA.txt',header=T,sep = "\t",row.names=1)
data_TCGA$sample=rownames(data_TCGA)

rownames(TCGA)=gsub('-01A','',rownames(TCGA))

TCGA=TCGA[intersect(rownames(TCGA),rownames(data_TCGA)),]
data_TCGA=data_TCGA[intersect(rownames(data_TCGA),rownames(TCGA)),]

aa=cbind(TCGA,data_TCGA)

df=as.data.frame(aa)
df$'Total.Mutation.Rate'=log(df$'Total.Mutation.Rate')
p <- ggplot(df, aes(x = Subtype.x, y = Total.Mutation.Rate, color = Subtype.x)) +
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
	labs(title = "PRJEB25780", x = "Group", y = "Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c('MP1','MP2'),c('MP2','MP3'),c('MP1','MP3')), method='wilcox.test',
                            label = "p.signif")

aa$Percent_Tumor='too Low'
aa$Percent_Tumor[which(aa$'Percent.Tumor.Cells' >= 20)]='Low'
aa$Percent_Tumor[which(aa$'Percent.Tumor.Cells' >= 40)]='Mid'
aa$Percent_Tumor[which(aa$'Percent.Tumor.Cells' >= 60)]='High'
aa$Percent_Tumor[which(aa$'Percent.Tumor.Cells' >= 80)]='too High'

aa$simplicity='Low'
aa$simplicity[which(aa$abs >= 0.2)]='Mid'
aa$simplicity[which(aa$abs >= 0.4)]='High'

ha = HeatmapAnnotation(Simplicity = as.character(aa$'simplicity'),
					   Percent_Tumor = as.character(aa$'Percent_Tumor'),
					   TCGA_subtype = as.character(aa$'TCGA.Subtype'),
                       NC_subtype = as.character(aa$'Subgroup'), 
                       Lauren_Class = as.character(aa$'Lauren.Class'),
                       Subtype=as.character(aa$'Subtype'), col = list(TCGA_subtype = c("CIN" = "#BA3630", "EBV" = "#9F70A6",'GS'="#427D39",'MSI'='#70174F'),# 设置surstat颜色
                                                                    Simplicity = c("Low" = "#BA3630", "Mid" = "#9F70A6",'High'="#427D39"),
																	NC_subtype = c("EP" = "#CF812A", "MP" = "#5A92D2"),# 设置gender 颜色
                                                                    Subtype= c("MP2" = '#509285', "MP1" = '#BA3630',"MP3" = '#A2855C'),# 设置stage颜色
                                                                   Lauren_Class = c("Diffuse" = "#9F70A6", "Intestinal" = "#427D39",'Mixed'="#282177",'Unidentified'='#70174F'),
																   Percent_Tumor = c("too Low"="#282177","High" = "#BA3630", "Mid" = "#9F70A6","Low"="#427D39","too High"="#70174F")))
Heatmap(t(aa[,18:20]),cluster_columns = F,cluster_rows = F,show_column_names = F,name = "heat",top_annotation = ha,col = c('#5A8FCA','#F2F2F0','#E31A1C'))

