bytlib load R-4.0.2
bytlib load gcc
R

library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)

set.seed(520)
dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/STAD_P_E_T_N.rds')
DefaultAssay(dat) <- "RNA"
TPM=as.data.frame(dat[["RNA"]]@data)

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
sc=as.data.frame(sc)
result=NULL
data_meta=as.matrix(sc)

Subtype=NULL
for(i in 1:nrow(data_meta)){
	ch_max=pmax(data_meta[i,1],data_meta[i,2],data_meta[i,3])
	if(ch_max == data_meta[i,1]){
		Subtype=c(Subtype,'MP_1')
	}	
	if(ch_max == data_meta[i,2]){
		Subtype=c(Subtype,'MP_2')
	}	
	if(ch_max == data_meta[i,3]){
		Subtype=c(Subtype,'MP_3')
	}	
}
data_meta=as.data.frame(data_meta)	
data_meta$Subtype=Subtype

dat@meta.data=cbind(dat@meta.data,data_meta)

DimPlot(dat,group.by="Subtype",reduction='umap',pt.size=0.1,label = TRUE)

test<-dat@meta.data 
test=test[which(test$seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)),]   
test$cluster<-factor(test$seurat_clusters,levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))

colors=c('0'='#F4F0EB','2'='#F4F0EB','3'='#F4F0EB','4'='#F4F0EB','5'='#F4F0EB',
		 '6'='#F4F0EB','9'='#F4F0EB','10'='#F4F0EB','11'='#F4F0EB','12'='#F4F0EB','13'='#F4F0EB',
		 '1'='#D1832C','7'='#D1832C','8'='#D1832C')
p1 =ggplot(test,aes(x=cluster,y=MP_1)) +
	geom_violin(aes(fill=seurat_clusters),alpha=0.7,outlier.alpha=0,width=0.7)+	
	geom_boxplot(width=0.2,alpha=0.7,outlier.alpha=0)+     
	labs(x = '')+
	scale_y_continuous(name = "Score")+
	ggtitle("MP1 Score")+
	scale_fill_manual(values=colors)+
	theme_bw()+
	#stat_compare_means(comparisons=Class,label="p.signif")+
	theme(axis.text.x=element_text(colour="black",size=14),	#设置x轴刻度标签的字体属性
		axis.text.y=element_text(size=12,face="plain"), 		#设置x轴刻度标签的字体属性
		axis.title.y=element_text(size = 15,face="plain"),	#设置y轴的标题的字体属性

		plot.title = element_text(size=12,face="bold",hjust = 0.5), #设置总标题的字体属性
		panel.grid.major = element_blank()
)
pdf('/public/workspace/liuqzh/gastric_cancer/figure2/GSE150290_MP_1_Score.pdf',width=12,height=5)
print(p1)
dev.off()

test$Group='Other'
test$Group[which(test$seurat_clusters %in% c(7,8))]='MP1'
test$Group<-factor(test$Group,levels = c('MP1','Other'))
colors_MP1=c('MP1'='#D1832C','Other'='#F4F0EB')
# 生成箱形图并添加散点图层
p <-ggplot(test, aes(x = Group, y = MP_1)) +
	geom_boxplot(aes(fill=Group),outlier.shape = NA) +  # 不显示离群值，且箱形图无填充
	scale_fill_manual(values=colors_MP1)+
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE150290 scRNA", x = "Group", y = "MP1 Score")  # 添加标签
# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP1", "Other")), method='wilcox.test',
							label = "p.signif", label.y = c(1))
pdf('/public/workspace/liuqzh/gastric_cancer/figure2/GSE150290_MP_1_Other.pdf',width=3,height=6)
print(p)
dev.off()

####
colors=c('0'='#F4F0EB','1'='#F4F0EB','3'='#F4F0EB','4'='#F4F0EB','6'='#F4F0EB',
		 '7'='#F4F0EB','8'='#F4F0EB','10'='#F4F0EB','11'='#F4F0EB','12'='#F4F0EB','13'='#F4F0EB',
		 '2'='#5B93D4','5'='#5B93D4','9'='#5B93D4')
p1 =ggplot(test,aes(x=cluster,y=MP_2)) +
	geom_violin(aes(fill=seurat_clusters),alpha=0.7,outlier.alpha=0,width=0.7)+	
	geom_boxplot(width=0.2,alpha=0.7,outlier.alpha=0)+     
	labs(x = '')+
	scale_y_continuous(name = "Score")+
	ggtitle("MP2 Score")+
	scale_fill_manual(values=colors)+
	theme_bw()+
	#stat_compare_means(comparisons=Class,label="p.signif")+
	theme(axis.text.x=element_text(colour="black",size=14),	#设置x轴刻度标签的字体属性
		axis.text.y=element_text(size=12,face="plain"), 		#设置x轴刻度标签的字体属性
		axis.title.y=element_text(size = 15,face="plain"),	#设置y轴的标题的字体属性

		plot.title = element_text(size=12,face="bold",hjust = 0.5), #设置总标题的字体属性
		panel.grid.major = element_blank()
)
pdf('/public/workspace/liuqzh/gastric_cancer/figure2/GSE150290_MP_2_Score.pdf',width=12,height=5)
print(p1)
dev.off()

test$Group='Other'
test$Group[which(test$seurat_clusters %in% c(5,9))]='MP2'
test$Group<-factor(test$Group,levels = c('MP2','Other'))
colors_MP2=c('MP2'='#5B93D4','Other'='#F4F0EB')
# 生成箱形图并添加散点图层
p <-ggplot(test, aes(x = Group, y = MP_2)) +
	geom_boxplot(aes(fill=Group),outlier.shape = NA) +  # 不显示离群值，且箱形图无填充
	scale_fill_manual(values=colors_MP2)+
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE150290 scRNA", x = "Group", y = "MP2 Score")  # 添加标签
# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP2", "Other")), method='wilcox.test',
							label = "p.signif", label.y = c(1.5))
pdf('/public/workspace/liuqzh/gastric_cancer/figure2/GSE150290_MP_2_Other.pdf',width=3,height=6)
print(p)
dev.off()
####
colors=c('0'='#F4F0EB','1'='#F4F0EB','2'='#F4F0EB','5'='#F4F0EB','6'='#F4F0EB',
		 '7'='#F4F0EB','8'='#F4F0EB','9'='#F4F0EB','10'='#F4F0EB','12'='#F4F0EB','13'='#F4F0EB'
		 '3'='#A97598','4'='#A97598','11'='#A97598')
p1 =ggplot(test,aes(x=cluster,y=MP_3)) +
	geom_violin(aes(fill=seurat_clusters),alpha=0.7,outlier.alpha=0,width=0.7)+	
	geom_boxplot(width=0.2,alpha=0.7,outlier.alpha=0)+     
	labs(x = '')+
	scale_y_continuous(name = "Score")+
	ggtitle("MP3 Score")+
	scale_fill_manual(values=colors)+
	theme_bw()+
	#stat_compare_means(comparisons=Class,label="p.signif")+
	theme(axis.text.x=element_text(colour="black",size=14),	#设置x轴刻度标签的字体属性
		axis.text.y=element_text(size=12,face="plain"), 		#设置x轴刻度标签的字体属性
		axis.title.y=element_text(size = 15,face="plain"),	#设置y轴的标题的字体属性

		plot.title = element_text(size=12,face="bold",hjust = 0.5), #设置总标题的字体属性
		panel.grid.major = element_blank()
)
pdf('/public/workspace/liuqzh/gastric_cancer/figure2/GSE150290_MP_3_Score.pdf',width=12,height=5)
print(p1)
dev.off()

test$Group='Other'
test$Group[which(test$seurat_clusters %in% c(3,4,11))]='MP3'
test$Group<-factor(test$Group,levels = c('MP3','Other'))
colors_MP3=c('MP3'='#A97598','Other'='#F4F0EB')
# 生成箱形图并添加散点图层
p <-ggplot(test, aes(x = Group, y = MP_3)) +
	geom_boxplot(aes(fill=Group),outlier.shape = NA) +  # 不显示离群值，且箱形图无填充
	scale_fill_manual(values=colors_MP3)+
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "GSE150290 scRNA", x = "Group", y = "MP3 Score")  # 添加标签
# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP3", "Other")), method='wilcox.test',
							label = "p.signif", label.y = c(1))
pdf('/public/workspace/liuqzh/gastric_cancer/figure2/GSE150290_MP_3_Other.pdf',width=3,height=6)
print(p)
dev.off()