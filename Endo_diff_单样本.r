#######
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(dplyr)
library(ggplot2)
library(ggtern)
library(GSVA)
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ComplexHeatmap)

set.seed(520)
#读取scRNA
scdata_score<-readRDS('/public/workspace/liuqzh/gastric_cancer/scalop-master/scdata_score.rds')

dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(dat) <- "RNA"

dat_cell=subset(dat,cells=which(dat$celltype %in% c('Endothelial')))
dat_cell@meta.data=dat_cell@meta.data[intersect(rownames(dat_cell@meta.data),rownames(scdata_score)),]
scdata_score=as.data.frame(scdata_score[intersect(rownames(dat_cell@meta.data),rownames(scdata_score)),])

dat_cell@meta.data$Endo_MP_1=scdata_score$Endo_MP_1
dat_cell@meta.data$Endo_MP_2=scdata_score$Endo_MP_2
dat_cell@meta.data$Endo_MP_3=scdata_score$Endo_MP_3
dat_cell@meta.data$Endo_MP_4=scdata_score$Endo_MP_4
dat_cell@meta.data$Endo_MP_5=scdata_score$Endo_MP_5
dat_cell@meta.data$Endo_MP_6=scdata_score$Endo_MP_6
dat_cell@meta.data$Endo_MP_7=scdata_score$Endo_MP_7

data_meta=as.data.frame(dat_cell@meta.data)
Subtype=NULL

for(i in 1:nrow(data_meta)){
	ch_max=max(data_meta[i,'Endo_MP_1'],data_meta[i,'Endo_MP_2'],data_meta[i,'Endo_MP_3'],
			   data_meta[i,'Endo_MP_4'],data_meta[i,'Endo_MP_5'],data_meta[i,'Endo_MP_6'],
			   data_meta[i,'Endo_MP_7'])
tmp_dat=ifelse(data_meta[i,'Endo_MP_1'] < ch_max,
			ifelse(data_meta[i,'Endo_MP_2'] < ch_max,
				ifelse(data_meta[i,'Endo_MP_3'] < ch_max,
					ifelse(data_meta[i,'Endo_MP_4'] < ch_max,
						ifelse(data_meta[i,'Endo_MP_5'] < ch_max,
							ifelse(data_meta[i,'Endo_MP_6'] < ch_max,
									'Endo_MP_7',
									'Endo_MP_6'),
									'Endo_MP_5'),
									'Endo_MP_4'),
									'Endo_MP_3'),
									'Endo_MP_2'),
									'Endo_MP_1')
Subtype=c(Subtype,tmp_dat)
}

Subtype=as.matrix(Subtype)
rownames(Subtype)=rownames(data_meta)
colnames(Subtype)='Subtype'

dat_cell@meta.data$Subtype=Subtype

pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_score_seurat.pdf',height=7,width=7)
DimPlot(dat_cell,group.by='Subtype')
dev.off()


list <- SplitObject(dat_cell, split.by = "Subtype")
sc.list=list()
sc.list = c(list$Endo_MP_1,list$Endo_MP_2,
			list$Endo_MP_3,list$Endo_MP_4,
			list$Endo_MP_5,list$Endo_MP_6,list$Endo_MP_7)

for(i in 1:length(sc.list)){
 tmp <- NormalizeData(sc.list[[i]])
 tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
 sc.list[[i]] <- tmp
}

anchors <- FindIntegrationAnchors(object.list = sc.list, dims = 1:30)
dat <- IntegrateData(anchorset = anchors, dims = 1:30)
dat <- ScaleData(dat, features = rownames(dat))
dat <- RunPCA(dat, features = VariableFeatures(object = dat))
dat <- FindNeighbors(dat, dims = 1:30)
dat <- FindClusters(dat, resolution = 1)
dat <- RunUMAP(dat, dims = 1:30)
dat <- RunTSNE(dat, dims = 1:30)
#saveRDS(dat, file="/public/workspace/liuqzh/gastric_cancer/Endo_score/dat_seurat_cluster.rds")
dat <- readRDS("/public/workspace/liuqzh/gastric_cancer/Endo_score/dat_seurat_cluster.rds")

DimPlot(dat,group.by="clas",reduction='umap',label=T)
DimPlot(dat,group.by="seurat_clusters",reduction='umap',label=T)
DimPlot(dat,group.by='Subtype')

Endo_clas=as.matrix(table(dat$clas,dat$Subtype))
STAD<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(STAD) <- "RNA"
MP_clas=as.matrix(table(STAD$clas,STAD$Subtype))

write.table(Endo_clas,'/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_clas.txt',row.names=T,col.names=T,quote=F,sep="\t")
write.table(MP_clas,'/public/workspace/liuqzh/gastric_cancer/Endo_score/MP_clas.txt',row.names=T,col.names=T,quote=F,sep="\t")

Endo_clas=read.table('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_clas.txt',header=T)
MP_clas=read.table('/public/workspace/liuqzh/gastric_cancer/Endo_score/MP_clas.txt',header=T)

Endo_clas$Sample=rownames(Endo_clas)
MP_clas$Sample=rownames(MP_clas)

Endo_MP_clas=merge(Endo_clas,MP_clas,'Sample')
rownames(Endo_MP_clas)=Endo_MP_clas$Sample
Endo_MP_clas=Endo_MP_clas[,-1]

data_f=Endo_MP_clas
p1=ggplot(data=data_f, aes(x=Endo_MP_5, y=MP_2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x = Endo_MP_5, y = MP_2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 0.5, label.y = 1) +  # 添加皮尔逊相关系数
	theme_minimal() +  # 使用简洁主题
	labs(
		title = "Relationship between Endo_MP_5 and MP_2",
		x = "Endo_MP_5",
		y = "MP_2"
	) +  # 添加图形标题和轴标题
	theme(
		plot.title = element_text(hjust = 0.5),  # 标题居中
		axis.text = element_text(color = "gray20"),  # 轴文字颜色
		axis.title = element_text(color = "gray20"),  # 轴标题颜色
		panel.grid.major = element_blank(),  # 移除主要网格线
		panel.grid.minor = element_blank(),  # 移除次要网格线
		panel.background = element_rect(fill = "white", color = "gray50")  # 背景颜色
	)

pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_5_MP_2.pdf',height=7,width=7)
p1
dev.off()

########
#saveRDS(dat_cell, "/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_score_data.rds")

###monocle
library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)

Endo_score_data <- readRDS("/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_score_data.rds")
DefaultAssay(Endo_score_data)<-'RNA'

STAD_Harmony_Epithelium=Endo_score_data
data=as.matrix(STAD_Harmony_Epithelium@assays$RNA@counts)
pd=new('AnnotatedDataFrame', data = STAD_Harmony_Epithelium@meta.data)
fData=data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd=new('AnnotatedDataFrame', data = fData)

dat <- newCellDataSet(data,phenoData = pd,featureData = fd,expressionFamily = negbinomial.size())

dat <- estimateSizeFactors(dat)
dat <- estimateDispersions(dat)
dat <- detectGenes(dat, min_expr = 1)
fData(dat)$use_for_ordering <- fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <-  row.names(subset(fData(dat), num_cells_expressed >= 10))

dat <- reduceDimension(dat, max_components = 2, num_dim = 20, verbose = T)
#dat <- clusterCells(dat, verbose = F, num_clusters = 2)
clustering_DEG_genes <- differentialGeneTest(dat[expressed_genes,], fullModelFormulaStr = '~Subtype', cores = 16)
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

#saveRDS(ordering_genes, "/public/workspace/liuqzh/gastric_cancer/Endo_score/mycds_Endo_score_data_ordering_genes.rds")
#ordering_genes<-readRDS('/public/workspace/liuqzh/gastric_cancer/Endo_score/mycds_Endo_score_data_ordering_genes.rds')

dat <- setOrderingFilter(dat, ordering_genes = ordering_genes)
dat <- reduceDimension(dat, method = 'DDRTree')

dat <- orderCells(dat)

mycds=dat
saveRDS(mycds, "/public/workspace/liuqzh/gastric_cancer/Endo_score/mycds_Endo_score_data_SeuratFeatures.rds")
mycds<-readRDS('/public/workspace/liuqzh/gastric_cancer/Endo_score/mycds_Endo_score_data_SeuratFeatures.rds')

colors=c("Endo_MP_1"="#BA3731","Endo_MP_2"="#A2855B","Endo_MP_3"="#5891D2",
		 "Endo_MP_4"="#CF812A","Endo_MP_5"="#509285","Endo_MP_6"="#427E39",
		 "Endo_MP_7"="#9F70A7")

##State轨迹分布图
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_score_data_monocle_subtype.pdf',width=6,height=6,useDingbats=F)
plot_cell_trajectory(mycds, color_by="Subtype",cell_size = 1)+scale_colour_manual(values=colors)
dev.off()

pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_score_data_monocle_Pseudotime.pdf',width=6,height=6,useDingbats=F)
plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 1)
dev.off()

###
Idents(dat_cell)=dat_cell$Subtype
Endo_marker=FindAllMarkers(object = dat_cell, only.pos = TRUE, logfc.threshold = 0.2)
gene_list=Endo_marker %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)

#write.table(gene_list,'/public/workspace/liuqzh/gastric_cancer/Endo_score/gene_list_Endo.txt',row.names=T,col.names=T,quote=F,sep="\t")

###############
###############
#######
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(dplyr)
library(ggplot2)
library(ggtern)
library(GSVA)
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ComplexHeatmap)

set.seed(520)
#读取scRNA
scdata_score<-readRDS('/public/workspace/liuqzh/gastric_cancer/scalop-master/scdata_score.rds')

dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(dat) <- "RNA"

dat_cell=subset(dat,cells=which(dat$celltype %in% c('Endothelial')))
dat_cell@meta.data=dat_cell@meta.data[intersect(rownames(dat_cell@meta.data),rownames(scdata_score)),]
scdata_score=as.data.frame(scdata_score[intersect(rownames(dat_cell@meta.data),rownames(scdata_score)),])

dat_cell@meta.data$Endo_MP_1=scdata_score$Endo_MP_1
dat_cell@meta.data$Endo_MP_2=scdata_score$Endo_MP_2
dat_cell@meta.data$Endo_MP_3=scdata_score$Endo_MP_3
dat_cell@meta.data$Endo_MP_4=scdata_score$Endo_MP_4
dat_cell@meta.data$Endo_MP_5=scdata_score$Endo_MP_5
dat_cell@meta.data$Endo_MP_6=scdata_score$Endo_MP_6
dat_cell@meta.data$Endo_MP_7=scdata_score$Endo_MP_7

data_meta=as.data.frame(dat_cell@meta.data)
Subtype=NULL

for(i in 1:nrow(data_meta)){
	ch_max=max(data_meta[i,'Endo_MP_1'],data_meta[i,'Endo_MP_2'],data_meta[i,'Endo_MP_3'],
			   data_meta[i,'Endo_MP_4'],data_meta[i,'Endo_MP_5'],data_meta[i,'Endo_MP_6'],
			   data_meta[i,'Endo_MP_7'])
tmp_dat=ifelse(data_meta[i,'Endo_MP_1'] < ch_max,
			ifelse(data_meta[i,'Endo_MP_2'] < ch_max,
				ifelse(data_meta[i,'Endo_MP_3'] < ch_max,
					ifelse(data_meta[i,'Endo_MP_4'] < ch_max,
						ifelse(data_meta[i,'Endo_MP_5'] < ch_max,
							ifelse(data_meta[i,'Endo_MP_6'] < ch_max,
									'Endo_MP_7',
									'Endo_MP_6'),
									'Endo_MP_5'),
									'Endo_MP_4'),
									'Endo_MP_3'),
									'Endo_MP_2'),
									'Endo_MP_1')
Subtype=c(Subtype,tmp_dat)
}

Subtype=as.matrix(Subtype)
rownames(Subtype)=rownames(data_meta)
colnames(Subtype)='Subtype'

dat_cell@meta.data$Subtype=Subtype
STAD<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(STAD) <- "RNA"

dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(dat) <- "RNA"
dat_all=subset(dat,cells=which(dat$celltype %in% c('Endothelial','Epithelium')))

dat_cell@meta.data$cell_names=colnames(dat_cell)
STAD@meta.data$cell_names=colnames(STAD)

cell_names=rbind(dat_cell@meta.data[c('Subtype','cell_names')],STAD@meta.data[c('Subtype','cell_names')])
cell_names=cell_names[intersect(rownames(dat_all@meta.data),rownames(cell_names)),]

dat_all@meta.data$Subtype=cell_names$Subtype

pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_PDGFRB.pdf',height=6,width=6)
FeaturePlot(object = dat_all, features = c('PDGFRB'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()

pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_OLFML2A.pdf',height=6,width=6)
FeaturePlot(object = dat_all, features = c('OLFML2A'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()

pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_BGN.pdf',height=6,width=6)
FeaturePlot(object = dat_all, features = c('BGN'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
#####单样本#####
dat_cell_single=subset(dat_all,cells=which(dat_all$clas %in% c('t1_5846')))
dat_cell_single=subset(dat_cell_single,cells=which(dat_cell_single$Subtype %in% c('MP_2','Endo_MP_1','Endo_MP_2','Endo_MP_3','Endo_MP_4','Endo_MP_5','Endo_MP_6','Endo_MP_7')))

dat_cell_Endo=subset(dat_cell_single,cells=which(dat_cell_single$Subtype %in% c('Endo_MP_5')))
expr_data <- as.data.frame(GetAssayData(dat_cell_Endo[["RNA"]], slot='data'))
#####
ncol(expr_data)
#72
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFRB'),]
length(which(dat_cell_Endo_data > 0))
#12
#####
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFA'),]
length(which(dat_cell_Endo_data > 0))
#23
#####
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFB'),]
which(dat_cell_Endo_data > 0)
length(which(dat_cell_Endo_data > 0))
#2

dat_cell_Endo=subset(dat_cell_single,cells=which(dat_cell_single$Subtype %in% c('MP_2')))
expr_data <- as.data.frame(GetAssayData(dat_cell_Endo[["RNA"]], slot='data'))

#####
ncol(expr_data)
#139
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFRB'),]
which(dat_cell_Endo_data > 0)
length(which(dat_cell_Endo_data > 0))
#0
#####
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFA'),]
which(dat_cell_Endo_data > 0)
length(which(dat_cell_Endo_data > 0))
#51
#####
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFB'),]
which(dat_cell_Endo_data > 0)
length(which(dat_cell_Endo_data > 0))
#0

colors=c("Endo_MP_1"="#BA3731","Endo_MP_2"="#A2855B","Endo_MP_3"="#5891D2",
		 "Endo_MP_4"="#CF812A","Endo_MP_5"="#509285","Endo_MP_6"="#427E39",
		 "Endo_MP_7"="#9F70A7","MP_1"="#21B573","MP_2"="#701750","MP_3"="#282178")

pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_t1_5846_PDGFRB.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('PDGFRB'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_t1_5846_PDGFA.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('PDGFA'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_t1_5846_PDGFB.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('PDGFB'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_t1_5846_CD44.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('CD44'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_t1_5846_LGALS9.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('LGALS9'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()


#####单样本#####
dat_cell_single=subset(dat_all,cells=which(dat_all$clas %in% c('Tu13')))
dat_cell_single=subset(dat_cell_single,cells=which(dat_cell_single$Subtype %in% c('MP_2','Endo_MP_1','Endo_MP_2','Endo_MP_3',
																				  'Endo_MP_4','Endo_MP_5','Endo_MP_6','Endo_MP_7')))
dat_cell_Endo=subset(dat_cell_single,cells=which(dat_cell_single$Subtype %in% c('Endo_MP_5')))
expr_data <- as.data.frame(GetAssayData(dat_cell_Endo[["RNA"]], slot='data'))
#####
ncol(expr_data)
#72
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFRB'),]
length(which(dat_cell_Endo_data > 0))
#35
#####
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFA'),]
length(which(dat_cell_Endo_data > 0))
#13
#####
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFB'),]
length(which(dat_cell_Endo_data > 0))
#15

dat_cell_Endo=subset(dat_cell_single,cells=which(dat_cell_single$Subtype %in% c('MP_2')))
expr_data <- as.data.frame(GetAssayData(dat_cell_Endo[["RNA"]], slot='data'))

#####
ncol(expr_data)
#455
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFRB'),]
length(which(dat_cell_Endo_data > 0))
#2
#####
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFA'),]
length(which(dat_cell_Endo_data > 0))
#31
#####
dat_cell_Endo_data=expr_data[which(rownames(expr_data)=='PDGFB'),]
length(which(dat_cell_Endo_data > 0))
#0

colors=c("Endo_MP_1"="#BA3731","Endo_MP_2"="#A2855B","Endo_MP_3"="#5891D2",
		 "Endo_MP_4"="#CF812A","Endo_MP_5"="#509285","Endo_MP_6"="#427E39",
		 "Endo_MP_7"="#9F70A7","MP_1"="#21B573","MP_2"="#701750","MP_3"="#282178")

pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_Tu13_PDGFRB.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('PDGFRB'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_Tu13_PDGFA.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('PDGFA'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_Tu13_PDGFB.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('PDGFB'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_Tu13_CD44.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('CD44'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()
pdf('/public/workspace/liuqzh/gastric_cancer/Endo_score/Endo_MP_Tu13_LGALS9.pdf',height=6,width=6)
FeaturePlot(object = dat_cell_single, features = c('LGALS9'),order=T,reduction='umap',cols = c("#a1a3a6","#BA3731"))
dev.off()


