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

Maligant_MP_1<-read.table('/public/workspace/liuqzh/gastric_cancer/Maligant_MP_1_TCGA.txt',header=T,row.names=1)
Maligant_MP_2<-read.table('/public/workspace/liuqzh/gastric_cancer/Maligant_MP_2_TCGA.txt',header=T,row.names=1)
Maligant_MP_3<-read.table('/public/workspace/liuqzh/gastric_cancer/Maligant_MP_3_TCGA.txt',header=T,row.names=1)
Maligant_MP_1$Rank_MP1=rank(Maligant_MP_1$Maligant_MP_1_norm)
Maligant_MP_2$Rank_MP2=rank(Maligant_MP_2$Maligant_MP_2_norm)
Maligant_MP_3$Rank_MP3=rank(Maligant_MP_3$Maligant_MP_3_norm)
TCGA_subtype=cbind(Maligant_MP_1,Maligant_MP_2,Maligant_MP_3)

data_meta=TCGA_subtype
Subtype=NULL
for(i in 1:nrow(data_meta)){
	ch_max=pmax(data_meta[i,4],data_meta[i,8],data_meta[i,12])
	TMP=ifelse(ch_max == data_meta[i,4],'MP_1',ifelse(ch_max == data_meta[i,8],'MP_2','MP_3'))
	Subtype=c(Subtype,TMP)
}
data_meta=as.data.frame(data_meta)	
data_meta$Subtype=Subtype

Maligant_MP_1=data_meta[which(data_meta$Subtype == 'MP_1'),]
Maligant_MP_2=data_meta[which(data_meta$Subtype == 'MP_2'),]
Maligant_MP_3=data_meta[which(data_meta$Subtype == 'MP_3'),]

#####MP1
ADNS=NULL
for(i in rownames(data_meta[which(data_meta$Subtype == 'MP_1'),])){
	Rank_MP2=data_meta[which(rownames(data_meta) == i),'Rank_MP2']
	Pval_MP2=data_meta[which(rownames(data_meta) == i),'Maligant_MP_2_pval']	
	MP2_tmp=data_meta[which(data_meta$Rank_MP2 < Rank_MP2),]
	ADNS_1=sum(MP2_tmp$Maligant_MP_2_pval)-(Rank_MP2-1)*Pval_MP2

	Rank_MP3=data_meta[which(rownames(data_meta) == i),'Rank_MP3']
	Pval_MP3=data_meta[which(rownames(data_meta) == i),'Maligant_MP_3_pval']	
	MP3_tmp=data_meta[which(data_meta$Rank_MP3 < Rank_MP3),]
	ADNS_2=sum(MP3_tmp$Maligant_MP_3_pval)-(Rank_MP3-1)*Pval_MP3
	
	ADNS_ALL=ADNS_1+ADNS_2
	ADNS=c(ADNS,ADNS_ALL)
}
Maligant_MP_1$ADNS=ADNS
#######
R0_MP_1=min(data_meta$Maligant_MP_1_pval)
Rn1_MP_1=max(data_meta$Maligant_MP_1_pval)
MP_1_ADDS=sum(data_meta$Maligant_MP_1_pval)+nrow(data_meta)*R0_MP_1
Maligant_MP_1=cbind(Maligant_MP_1,MP_1_ADDS)
Maligant_MP_1$Simplicity=(Maligant_MP_1$MP_1_ADDS-Maligant_MP_1$ADNS)*(Rn1_MP_1-R0_MP_1)/nrow(Maligant_MP_1)
write.table(Maligant_MP_1, "/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_1.txt", row.names = T, sep = "\t")


#####MP2
ADNS=NULL
for(i in rownames(data_meta[which(data_meta$Subtype == 'MP_2'),])){
	Rank_MP1=data_meta[which(rownames(data_meta) == i),'Rank_MP1']
	Pval_MP1=data_meta[which(rownames(data_meta) == i),'Maligant_MP_1_pval']	
	MP1_tmp=data_meta[which(data_meta$Rank_MP1 < Rank_MP1),]
	ADNS_1=sum(MP1_tmp$Maligant_MP_1_pval)-(Rank_MP1-1)*Pval_MP1

	Rank_MP3=data_meta[which(rownames(data_meta) == i),'Rank_MP3']
	Pval_MP3=data_meta[which(rownames(data_meta) == i),'Maligant_MP_3_pval']	
	MP3_tmp=data_meta[which(data_meta$Rank_MP3 < Rank_MP3),]
	ADNS_2=sum(MP3_tmp$Maligant_MP_3_pval)-(Rank_MP3-1)*Pval_MP3
	
	ADNS_ALL=ADNS_1+ADNS_2
	ADNS=c(ADNS,ADNS_ALL)
}
Maligant_MP_2$ADNS=ADNS
#######
R0_MP_2=min(data_meta$Maligant_MP_2_pval)
Rn1_MP_2=max(data_meta$Maligant_MP_2_pval)
MP_2_ADDS=sum(data_meta$Maligant_MP_2_pval)+nrow(data_meta)*R0_MP_2
Maligant_MP_2=cbind(Maligant_MP_2,MP_2_ADDS)
Maligant_MP_2$Simplicity=(Maligant_MP_2$MP_2_ADDS-Maligant_MP_2$ADNS)*(Rn1_MP_2-R0_MP_2)/nrow(Maligant_MP_2)
write.table(Maligant_MP_2, "/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_2.txt", row.names = T, sep = "\t")


#####MP3
ADNS=NULL
for(i in rownames(data_meta[which(data_meta$Subtype == 'MP_3'),])){
	Rank_MP1=data_meta[which(rownames(data_meta) == i),'Rank_MP1']
	Pval_MP1=data_meta[which(rownames(data_meta) == i),'Maligant_MP_1_pval']	
	MP1_tmp=data_meta[which(data_meta$Rank_MP1 < Rank_MP1),]
	ADNS_1=sum(MP1_tmp$Maligant_MP_1_pval)-(Rank_MP1-1)*Pval_MP1

	Rank_MP2=data_meta[which(rownames(data_meta) == i),'Rank_MP2']
	Pval_MP2=data_meta[which(rownames(data_meta) == i),'Maligant_MP_2_pval']	
	MP2_tmp=data_meta[which(data_meta$Rank_MP2 < Rank_MP2),]
	ADNS_2=sum(MP2_tmp$Maligant_MP_2_pval)-(Rank_MP2-1)*Pval_MP2
	
	ADNS_ALL=ADNS_1+ADNS_2
	ADNS=c(ADNS,ADNS_ALL)
}
Maligant_MP_3$ADNS=ADNS
#######
R0_MP_3=min(data_meta$Maligant_MP_3_pval)
Rn1_MP_3=max(data_meta$Maligant_MP_3_pval)
MP_3_ADDS=sum(data_meta$Maligant_MP_3_pval)+nrow(data_meta)*R0_MP_3
Maligant_MP_3=cbind(Maligant_MP_3,MP_3_ADDS)
Maligant_MP_3$Simplicity=(Maligant_MP_3$MP_3_ADDS-Maligant_MP_3$ADNS)*(Rn1_MP_3-R0_MP_3)/nrow(Maligant_MP_3)
write.table(Maligant_MP_3, "/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_3.txt", row.names = T, sep = "\t")

#################
normalize <- function(x) {
return((x - min(x)) / (max(x) - min(x)))
}

Maligant_MP_1<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_1.txt',header=T,row.names=1)
Maligant_MP_2<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_2.txt',header=T,row.names=1)
Maligant_MP_3<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_3.txt',header=T,row.names=1)
Maligant_MP_1=arrange(Maligant_MP_1,desc(Simplicity))
Maligant_MP_1$Simplicity_score=normalize(Maligant_MP_1$Simplicity)

Maligant_MP_2=arrange(Maligant_MP_2,desc(Simplicity))
Maligant_MP_2$Simplicity_score=normalize(Maligant_MP_2$Simplicity)

Maligant_MP_3=arrange(Maligant_MP_3,desc(Simplicity))
Maligant_MP_3$Simplicity_score=normalize(Maligant_MP_3$Simplicity)

colnames(Maligant_MP_1)[15]='ADDS'
colnames(Maligant_MP_2)[15]='ADDS'
colnames(Maligant_MP_3)[15]='ADDS'

DATA_MP=rbind(Maligant_MP_1,Maligant_MP_2,Maligant_MP_3)
DATA_MP$MP1 = -log10(DATA_MP$Maligant_MP_1_pval)
DATA_MP$MP2 = -log10(DATA_MP$Maligant_MP_2_pval)
DATA_MP$MP3 = -log10(DATA_MP$Maligant_MP_3_pval)

pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.pdf",width=12,height=3)
library(circlize)
DATA_MP$simplicity='Low'
DATA_MP$simplicity[which(DATA_MP$Simplicity_score >= 0.8)]='High'

ha = HeatmapAnnotation(Subtype=as.character(DATA_MP$'Subtype'),
					   simplicity = as.character(DATA_MP$'simplicity'),
					   col = list(simplicity = c("Low" = "#427D39",'High'="#BA3630"),
								  Subtype= c("MP_1"='#D1832C',"MP_2"='#5B92D3',"MP_3"='#A97597')))

mycols <- colorRamp2(breaks = c(0,1.5,5),
                    colors = c('#5A8FCA','#F2F2F0','#E31A1C'))
Heatmap(t(DATA_MP[,c(18:20)]),cluster_columns = F,cluster_rows = F,show_column_names = F,name = "heat",top_annotation = ha,
		col = mycols)
dev.off()

#############
data_TCGA=read.table('/public/workspace/liuqzh/gastric_cancer/info_data/TCGA.txt',header=T,sep = "\t",row.names=1)
data_TCGA$sample=rownames(data_TCGA)

TCGA=DATA_MP
rownames(TCGA)=gsub('-01A','',rownames(TCGA))

TCGA=TCGA[intersect(rownames(TCGA),rownames(data_TCGA)),]
data_TCGA=data_TCGA[intersect(rownames(TCGA),rownames(data_TCGA)),]

aa=cbind(TCGA,data_TCGA)

###
pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_all_result.pdf",width=12,height=3)
library(circlize)
df=aa
df$simplicity='Low'
df$simplicity[which(df$Simplicity_score >= 0.8)]='High'

df$TNM='none'
df$TNM[which(df$TNM.Stage %in% c('Stage_IA','Stage_IB'))]='Stage_I'
df$TNM[which(df$TNM.Stage %in% c('Stage_IIA','Stage_IIB'))]='Stage_II'
df$TNM[which(df$TNM.Stage %in% c('Stage_IIIA','Stage_IIIB','Stage_IIIC'))]='Stage_III'
df$TNM[which(df$TNM.Stage %in% c('Stage_IV'))]='Stage_IV'

df$T_Stage='none'
df$T_Stage[which(df$'Pathologic.T' %in% c('T1a','T1b'))]='T1'
df$T_Stage[which(df$'Pathologic.T' %in% c('T2'))]='T2'
df$T_Stage[which(df$'Pathologic.T' %in% c('T3'))]='T3'
df$T_Stage[which(df$'Pathologic.T' %in% c('T4','T4a','T4b'))]='T4'

df$N_Stage='none'
df$N_Stage[which(df$'Pathologic.N' %in% c('N0'))]='N0'
df$N_Stage[which(df$'Pathologic.N' %in% c('N1'))]='N1'
df$N_Stage[which(df$'Pathologic.N' %in% c('N2'))]='N2'
df$N_Stage[which(df$'Pathologic.N' %in% c('N3'))]='N3'

df$M_Stage='none'
df$M_Stage[which(df$'Pathologic.M' %in% c('M0'))]='M0'
df$M_Stage[which(df$'Pathologic.M' %in% c('M1'))]='M1'

ha = HeatmapAnnotation(Subtype=as.character(df$'Subtype'),
					   simplicity = as.character(df$'simplicity'),
					   TNM = as.character(df$'TNM'),
					   T_Stage = as.character(df$'T_Stage'),
					   N_Stage = as.character(df$'N_Stage'),
					   M_Stage = as.character(df$'M_Stage'),
					   col = list(simplicity = c("Low" = "#427D39",'High'="#BA3630"),
								  Subtype= c("MP_1"='#D1832C',"MP_2"='#5B92D3',"MP_3"='#A97597'),
								  TNM = c("Stage_I"='#D1832C',"Stage_II"='#5B92D3',"Stage_III"='#A97597',"Stage_IV"='#427D39',"none"='grey'),
								  T_Stage = c("T1"='#D1832C',"T2"='#5B92D3',"T3"='#A97597',"T4" = '#427D39',"none"='grey'),
								  N_Stage = c("N0"='#D1832C',"N1"='#5B92D3',"N2"='#A97597',"N3" = '#427D39',"none"='grey'),
								  M_Stage = c("M0"='#D1832C',"M1"='#5B92D3',"none"='grey')))

mycols <- colorRamp2(breaks = c(0,1.5,5),
                    colors = c('#5A8FCA','#F2F2F0','#E31A1C'))
Heatmap(t(df[,c(18:20)]),cluster_columns = F,cluster_rows = F,show_column_names = F,name = "heat",top_annotation = ha,
		col = mycols)
dev.off()

df=as.data.frame(aa)
colnames(df)[41]='old'
df$'Total.Mutation.Rate'=log(df$'Total.Mutation.Rate')
p <- ggplot(df, aes(x = Subtype, y = Total.Mutation.Rate, color = Subtype)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("High" = "red", "Low" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "TCGA-STAD", x = "Group", y = "Total.Mutation.Rate")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c('MP_1','MP_2'),c('MP_2','MP_3'),c('MP_1','MP_3')), method='wilcox.test',
                            label = "p.signif")

pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_Mutation_Rate.pdf",width=4,height=6)
print(p)
dev.off()

#####
df$Percent_Tumor='Low'
df$Percent_Tumor[which(df$'Percent.Tumor.Cells' >= 40)]='Mid'
df$Percent_Tumor[which(df$'Percent.Tumor.Cells' >= 60)]='High'

df$simplicity='Low'
df$simplicity[which(df$Simplicity_score >= 0.8)]='High'

ha = HeatmapAnnotation(Simplicity = as.character(df$'simplicity'),
					   Percent_Tumor = as.character(df$'Percent_Tumor'),
					   TCGA_subtype = as.character(df$'TCGA.Subtype'),
                       NC_subtype = as.character(df$'Subgroup'), 
                       Lauren_Class = as.character(df$'Lauren.Class'),
                       Subtype=as.character(df$'Subtype'), col = list(TCGA_subtype = c("CIN" = "#BA3630", "EBV" = "#9F70A6",'GS'="#427D39",'MSI'='#70174F'),# 设置surstat颜色
                                                                    Simplicity = c("Low" = "#427D39",'High'="#BA3630"),
																	NC_subtype = c("EP" = "#CF812A", "MP" = "#5A92D2"),# 设置gender 颜色
                                                                    Subtype= c("MP_2" = '#509285', "MP_1" = '#BA3630',"MP_3" = '#A2855C'),# 设置stage颜色
                                                                   Lauren_Class = c("Diffuse" = "#9F70A6", "Intestinal" = "#427D39",'Mixed'="#282177",'Unidentified'='#70174F'),
																   Percent_Tumor = c("High" = "#BA3630", "Mid" = "#9F70A6","Low"="#5A92D2")))

pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_subtype.pdf",width=12,height=3)
Heatmap(t(df[,18:20]),cluster_columns = F,cluster_rows = F,show_column_names = F,name = "heat",top_annotation = ha,col = c('#5A8FCA','#F2F2F0','#E31A1C'))
dev.off()


