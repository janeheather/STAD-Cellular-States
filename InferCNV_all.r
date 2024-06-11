bytlib load R-4.0.2
bytlib load gcc
bytlib load JAGS-4.3.0
R

library(readr)
library(infercnv)
library(stringr)
library(dplyr)
library(circlize)
library(tidyverse)
library(ComplexHeatmap)
library(Seurat)
library(RColorBrewer)
full_join1 <- function(x,y) {full_join(x, y, by = c("id"))}

samples <- c('t1_5846','t1_5866','t1_5931','t1_6207','t1_6342','t1_6592','t2_5866','t2_5931',
			 'CRA002586-D01-T','CRA002586-D02-T','CRA002586-D03-T','CRA002586-P01-T','CRA002586-P03-T','CRA002586-P04-T','CRA002586-P05-T','CRA002586-P06-T',
			 'P01_EGC','P02_AGC','P04_EGC','P05_EGC','P06_EGC','P07_EGC','P08_EGC','P09_EGC',
			 'P10_EGC','P11_AGC','P12_AGC','P13_EGC','P15_AGC','P16_EGC','P17_EGC','P18_EGC',
			 'P19_EGC','P20_AGC','P21_EGC','P23_AGC',
			 'GSE167297-P1-1','GSE167297-P1-2','GSE167297-P2-1','GSE167297-P3-1',
			 'Tu1','Tu10','Tu11','Tu12','Tu13','Tu15','Tu16','Tu17','Tu18','Tu19','Tu2','Tu20','Tu21',
			 'Tu22','Tu23','Tu24','Tu25','Tu26','Tu3','Tu4','Tu5','Tu6','Tu8','Tu9')

expr.datas <- sapply(samples, function(x) {
	infercnv_obj = readRDS(as.character(str_glue("/public/workspace/liuqzh/gastric_cancer/Infercnv_all/", x, "/run.final.infercnv_obj")))
	tmp <- as.data.frame(infercnv_obj@expr.data[,infercnv_obj@observation_grouped_cell_indices$'Epithelium'])
	#rownames(tmp) <- rownames(dat@expr.data)
	tmp$id <- rownames(tmp)
	tmp
})
res <- Reduce(full_join1, expr.datas)
rownames(res) <- res$id
res$id <- NULL
res[is.na(res)] <- 1

geneInfo <- read.table("/public/workspace/liuqzh/lqz_ref/GRCh38_position.txt", sep = "\t", header = F)
colnames(geneInfo) <- c("Gene", "Chr", "Start", "End")
geneInfo$Chr <- as.character(geneInfo$Chr)
geneInfo <- geneInfo[geneInfo$Chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", 
										 "9", "10", "11", "12", "13", "14", "15",
										 "16", "17", "18", "19", "20", "21", "22"), ]
geneInfo$Chr <- factor(geneInfo$Chr,levels=c("1", "2", "3", "4", "5", "6", "7", "8", 
											 "9", "10", "11", "12", "13", "14", "15",
											 "16", "17", "18", "19", "20", "21", "22")) 
rownames(geneInfo) <- geneInfo$Gene
geneInfo1 <- geneInfo[intersect(geneInfo$Gene, rownames(res)),]
geneInfo1 <- geneInfo1[!is.na(geneInfo1$Chr), ]
geneInfo1 <- geneInfo1 %>% dplyr::arrange(Chr, Start, End)
geneInfo1 <- as.data.frame(geneInfo1)
rownames(geneInfo1) <- geneInfo1$Gene
geneInfo2 <- as.data.frame(geneInfo1[,2])
colnames(geneInfo2) <- "Chr"

cellInfo <- data.frame(
	Sample = as.character(sapply(colnames(res), function(x) {
		if (str_count(x, "_") == 2) {
			str_glue(strsplit(x, "_")[[1]][1], "_", strsplit(x, "_")[[1]][2])
		} else {
			strsplit(x, "_")[[1]][1]
		}
		}
	)),
	row.names = colnames(res)
)

cellInfo_II=cellInfo[which(grepl('^Tu',cellInfo$Sample)),,drop=F]
strsplit(cellInfo_II$Sample, "_")[[1]][1]

cellInfo_II$Sample=sub('_.*','',cellInfo_II$Sample)
cellInfo_I=cellInfo[-which(grepl('^Tu',cellInfo$Sample)),,drop=F]

cellInfo=rbind(cellInfo_I,cellInfo_II)
################

res1 <- res[as.character(geneInfo1$Gene),]
geneInfo3 <- as.character(geneInfo2$Chr)

left_annotation <- rowAnnotation(df = cellInfo)
top_annotation <- HeatmapAnnotation(df = geneInfo3)

pdf("/public/workspace/liuqzh/gastric_cancer/landscape_infercnv.pdf", width = 15,height = 18, useDingbats = F)
ht <- ComplexHeatmap::Heatmap(
	t(res1), 
	top_annotation = top_annotation,
	left_annotation = left_annotation,
	name = "inferCNV", 
	col = colorRamp2(c(0.9, 1, 1.1), c('blue', "white", 'red')),
	row_title_rot = 0, column_title_rot = 90,
	show_row_names = FALSE, show_column_names = FALSE, 
	cluster_rows = FALSE, cluster_columns = FALSE,
	row_split = cellInfo$Sample,
	column_split = factor(geneInfo3, levels = unique(geneInfo3)), 
	#column_gap = unit(0, "mm"),
	column_title_gp = gpar(fontsize = 6),
	#column_labels = clabs,
	column_names_gp = gpar(fontsize = 6),
	border = TRUE
)
draw(ht)
dev.off()

#########

inferCNV_data=as.data.frame((t(res1)-1)^2)
ID=unique(cellInfo$Sample)
Chr=unique(geneInfo3)

result=NULL
for(i in ID){
	print(i)
	
	tmp_Sample_CNV=inferCNV_data[which(rownames(inferCNV_data) %in% rownames(cellInfo[which(cellInfo$Sample==i),,drop=F])),]
	Sample_CNV_mean_cell=nrow(tmp_Sample_CNV)
	
	Sample_CNV_chr=NULL
	for(j in Chr){
		tmp_Sample_CNV_chr=tmp_Sample_CNV[,which(geneInfo3==j)]
		tmp_Chr_gene_avg=as.matrix(apply(tmp_Sample_CNV_chr,2,sum))/Sample_CNV_mean_cell
		Chr_score=as.matrix(sum(tmp_Chr_gene_avg)/length(tmp_Chr_gene_avg))
		rownames(Chr_score)=paste('Chr_',j,sep='')
		Sample_CNV_chr=rbind(Sample_CNV_chr,Chr_score)
	}
	colnames(Sample_CNV_chr)=i	
	result=cbind(result,Sample_CNV_chr)
}

library(ggplot2)
library(ggpubr)
library(reshape2)

df_long <- melt(result, id.vars = colnames(result), measure.vars = rownames(result), variable.name = "group", value.name = "Value")
df=df_long
colnames(df)[1]='Chr'
colnames(df)[2]='Sample'
df$group <- factor(df$Chr, levels = c("Chr_1","Chr_2","Chr_3","Chr_4","Chr_5",
									  "Chr_6","Chr_7","Chr_8","Chr_9","Chr_10",
									  "Chr_11","Chr_12","Chr_13","Chr_14","Chr_15",
									  "Chr_16","Chr_17","Chr_18","Chr_19","Chr_20",
									  "Chr_21","Chr_22"))  # 设置分组的顺序


# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = Value, color = group)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "ScRNA - seq inferCNV score", x = "Group", y = "Infercnv Score")  # 添加标签

pdf('/public/workspace/liuqzh/gastric_cancer/InferCNV_Score/InferCNV_Score.pdf',height=6,width=12)
print(p)
dev.off()	