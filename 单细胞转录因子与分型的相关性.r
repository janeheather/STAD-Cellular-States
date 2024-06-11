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

Cho_Epi=read.csv("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/Cho_Epithelium.auc_mtx.csv",sep=",",header=T,row.names=1)
colnames(Cho_Epi)=gsub("\\...","",colnames(Cho_Epi))

STAD_P<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(STAD_P)<-'RNA'
dat_ms=subset(STAD_P,cells=which(STAD_P$celltype=='Epithelium'))
DefaultAssay(dat_ms) <- "RNA"

Module_Sample_Subtype_III<-readRDS("/public/workspace/liuqzh/gastric_cancer/Module_Sample_Subtype_III.rds")
data_meta=Module_Sample_Subtype_III

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
Subtype=as.matrix(Subtype)
rownames(Subtype)=rownames(Module_Sample_Subtype_III)
colnames(Subtype)='Subtype'
dat_ms@meta.data$Subtype=Subtype

Cho_Epi=Cho_Epi[intersect(rownames(Cho_Epi),rownames(dat_ms@meta.data)),]
Cho_Epi=t(Cho_Epi)

data=dat_ms
data[['tfswmean']] <- Cho_Epi %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "tfswmean"

# Scale the data
data <- ScaleData(data)
data@assays$tfswmean@data <- data@assays$tfswmean@scale.data
dat1=data

rssMat <- sapply(rownames(dat1), function(x) {
	sapply(c("MP_1","MP_2","MP_3"), function(y) {
		tmp = rbind(
			value = as.numeric(dat1@assays$tfswmean@data[x,] + abs(min(dat1@assays$tfswmean@data))), 
			groups = as.numeric(as.character(ifelse(as.character(dat1$'Subtype') == y, 1, 0)))
		)
		1 - JSD(tmp, unit = 'log2', est.prob = "empirical")
	})
})
rownames(rssMat) <- c("MP_1","MP_2","MP_3")
rssMat <- as.data.frame(t(rssMat))

write.csv(rssMat, "/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/cho_epi_bbknn_rssMat.csv")

plotlist <- lapply(c("MP_1","MP_2","MP_3"), function(x) {
	library(dplyr)
	library(ggrepel)
	library(ggthemes)
	library(gridExtra)
	dat <- data.frame(
		rss = as.numeric(rssMat[[x]]),
		tfs = rownames(rssMat)
	)
	dat1 <- dat %>% arrange(desc(rss))
	dat1$index <- 1:nrow(dat1)
	dat1$showtfs <- ifelse(dat1$index <= 20, "1", "0")
	dat2 <- dat1[1:20,]
	library(ggplot2)
	#library(ggprism)

	p <- ggplot(data = dat1[1:50,], aes(x = index, y = rss, label=tfs)) + 
		geom_point(aes(color = showtfs)) + theme_classic() + 
		scale_color_manual(values = c("grey", "red")) +
		geom_point(data = dat2, alpha = 1, size = 2, shape = 1, stroke = 1, color = "black") +
		geom_text_repel(
			data = dat2, show.legend = FALSE, #不显示图例
			size = 8, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")
		) +
		guides(color=guide_legend(title = NULL)) +
		#theme_prism(base_size=12) + 
		xlab("tfs") + ylab("RSS") + ggtitle(x) +
		theme(legend.position = "none") +
		theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())
	return(p)
})

pdf("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/cho_epi_rssMat.pdf", useDingbats = F, height = 5, width = 5)
cowplot::plot_grid(plotlist = plotlist[1],ncol=1, nrow = 1)
cowplot::plot_grid(plotlist = plotlist[2],ncol=1, nrow = 1)
cowplot::plot_grid(plotlist = plotlist[3],ncol=1, nrow = 1)
dev.off()

Cho_Epi=t(Cho_Epi)
M=cor(Cho_Epi,method='pearson')
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0", "#FFFFFF",  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B","#67001F" ))

pdf("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/pySCENIC_heatmap.pdf", useDingbats = F, height = 6, width = 6)
p1 = corrplot(M, method = 'color',addrect = 3,
			order = 'hclust', 
			col=col2(100), diag = FALSE,sig.level = 0.05,
			hclust.method="complete",cl.pos = "n", tl.pos = "n")
dev.off()

#saveRDS(M,"/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/pySCENIC_M.rds")

###########转录因子活性的相关性热图#############
data_cor=cbind(data_meta,Cho_Epi)
data_cor=cor(data_cor,method='pearson')
data_cor=data_cor[,-c(1:3)]

data_cor=data_cor[,intersect(colnames(p$corr),colnames(data_cor))]
data_cor_plot=data_cor[c(1:3),]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0", "#FFFFFF",  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B","#67001F" ))

pdf("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/pySCENIC_MP1_MP2_MP3_cor.pdf", useDingbats = F, height = 3, width = 6)
pheatmap(data_cor_plot,cluster_rows = FALSE,cluster_cols = FALSE,col=col2(100),show_rownames=FALSE,show_colnames=FALSE)
dev.off()

