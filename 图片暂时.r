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
library(ggplot2)
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

STAD_FPKM=readRDS("/public/workspace/liuqzh/gastric/STAD_FPKM_data_exp.rds")
STAD_FPKM=STAD_FPKM[,grep('01A',colnames(STAD_FPKM))]
colnames(STAD_FPKM)=gsub('-01A','',colnames(STAD_FPKM))

TPM=as.data.frame(STAD_FPKM)
TPM<-TPM[rowSums(TPM)>0,]
TPM<-as.matrix(TPM)

genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/subtype_III.csv",sep = ",",row.names = 1,header = F)
geneset<-c()
for (i in 1:nrow(genelist)) {
	aa<-genelist[i,]
	aa<-aa[!is.na(aa)]
	aa<-aa[!duplicated(aa)]
	aa<-aa[aa%in%rownames(TPM)]
	aa<-list(as.character(as.matrix(aa)))
	geneset<-c(geneset,aa)
}
names(geneset)<-rownames(genelist)

colCenter = function(m, by = 'mean') {
	m = as.matrix(m)
	if (by == 'mean')  by = T
	else if (by == 'median') by = matrixStats::colMedians(m)
	else stopifnot(is.numeric(by) & length(by) == ncol(m))
	scale(m, center = by, scale = F)
}
sc = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

data_meta=sc
Subtype=NULL
for(i in 1:nrow(data_meta)){
	ch_max=pmax(data_meta[i,1],data_meta[i,2],data_meta[i,3])
	if(ch_max == data_meta[i,1]){
		Subtype=c(Subtype,'MP1')
	}
	if(ch_max == data_meta[i,2]){
		Subtype=c(Subtype,'MP2')
	}	
	if(ch_max == data_meta[i,3]){
		Subtype=c(Subtype,'MP3')
	}	
}
data_meta=as.data.frame(data_meta)	
data_meta$Subtype=Subtype

library(ComplexHeatmap)
Class_inf<-read.csv('/public/workspace/liuqzh/gastric_cancer/TCGA_subtype/subype_III_and_subtype_classical_order.csv',header=T,row.names=1)
Subtype_III_tcga=data_meta

TIDE_TCGA_STAD<-read.csv('/public/workspace/liuqzh/gastric_cancer/TIDE/TCGA_STAD_RNASeq_norm_subtract.txt',header=T,row.names=1,sep='\t')

Class_inf=Class_inf[,8:10]
Class_inf$Sample_ID=rownames(Class_inf)
Subtype_III_tcga$Sample_ID=rownames(Subtype_III_tcga)
a=merge(Class_inf,Subtype_III_tcga,'Sample_ID')

rownames(a)=a$Sample_ID
a=a[,-1]
a=a[which(a$Subtype != 'Other'),]

MP1=a[which(a$Subtype == 'MP1'),]
MP1=dplyr::arrange(MP1,desc(MP_1))
MP2=a[which(a$Subtype == 'MP2'),]
MP2=dplyr::arrange(MP2,desc(MP_2))
MP3=a[which(a$Subtype == 'MP3'),]
MP3=dplyr::arrange(MP3,desc(MP_3))

a=rbind(MP1,MP2,MP3)

ha = HeatmapAnnotation(TCGA_subtype = as.character(a$TCGA_subtype),
                       NC_subtype = as.character(a$NC_subtype), 
                       Lauren_Class = as.character(a$Lauren_Class),
                       Subtype=as.character(a$Subtype), col = list(TCGA_subtype = c("CIN" = "red", "EBV" = "green",'GS'="black",'MSI'='pink'),# 设置surstat颜色
                                                                    NC_subtype = c("EP" = "yellow", "MP" = "blue"),# 设置gender 颜色
                                                                    Subtype= c("MP2" = 'grey', "MP1" = 'red',"MP3" = 'pink'),# 设置stage颜色
                                                                   Lauren_Class=c("Diffuse" = "red", "Intestinal" = "green",'Mixed'="black",'Unidentified'='pink')))


Heatmap(t(a[,c(4,5,6)]),cluster_columns = F,show_column_names = F,name = "heat",top_annotation = ha,col = c('#5A8FCA','#F2F2F0','#E31A1C'))

TIDE_TCGA_STAD$Sample=rownames(TIDE_TCGA_STAD)
data_meta$Sample=rownames(data_meta)

data_TIDE=merge(TIDE_TCGA_STAD,data_meta,'Sample')

data_TIDE$Subtype <- factor(data_TIDE$Subtype, levels = c("MP1", "MP2","MP3"))  #设置分组的顺序
rownames(data_TIDE)=data_TIDE[,1]
data_TIDE=data_TIDE[,-1]

df_long <- melt(data_TIDE, id.vars = c("Subtype"), measure.vars = c("MDSC", "CAF", "M2", "Exclusion", "Dysfunction"), variable.name = "group", value.name = "Value")

df=df_long
df$group <- factor(df$group, levels = c("MDSC", "CAF", "M2", "Exclusion", "Dysfunction"))  # 设置分组的顺序
df$Subtype <- factor(df$Subtype, levels = c("MP1","MP2","MP3"))  # 设置分组的顺序
df$Value = as.numeric(df$Value)

p <-ggplot(df, aes(x = group, y = Value, fill = Subtype)) +
	geom_boxplot(outlier.shape = NA) +  # 不显示离群值，且箱形图无填充
	scale_color_manual()+
	scale_fill_manual(values=c("MP1"="#D1832C","MP2"="#5B93D4","MP3"="#A97598")) +  # 设置颜色
	#geom_point(fill = Subtype,size=2,color='black') +
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "STAD_TIDE", x = "Group", y = "TIDE Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP1", "MP2")), method='wilcox.test',label = "p.signif")

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE206329/STAD_OVm_Data_MP2.pdf',height=6,width=3)
print(p)
dev.off()

#################
####MDSC
df_tmp=df[which(df$group=='MDSC'),]

p <-ggplot(df_tmp, aes(x = group, y = Value, fill = Subtype)) +
	geom_boxplot(outlier.shape = NA) +  # 不显示离群值，且箱形图无填充
	scale_color_manual()+
	scale_fill_manual(values=c("MP1"="#D1832C","MP2"="#5B93D4","MP3"="#A97598")) +  # 设置颜色
	#geom_point(fill = Subtype,size=2,color='black') +
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "STAD_TIDE", x = "Group", y = "TIDE Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP1", "MP2"),c("MP1", "MP3"),c("MP2", "MP3")),
							method='wilcox.test',label = "p.signif", label.y = c(2, 2,2.2))








