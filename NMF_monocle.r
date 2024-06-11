#############figure2_A#################
{
bytlib load R-4.0.2
bytlib load gcc
R

library(Seurat)
library(dplyr)
library(ggplot2)
DAT <- readRDS("/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/STAD_P_E_T_N.rds")
DefaultAssay(DAT)<-'RNA'

a<-table(DAT$orig.ident) %>% as.data.frame()
b<-a[,1]
c<-data.frame()
for(i in b){
	tmp<-subset(DAT,subset=orig.ident==i)
	tmp1<-table(tmp$subtype) %>% as.data.frame()
	tmp2<-cbind(tmp1,i)
	c<-rbind(c,tmp2)
}
colnames(c)<-c("subtype","percent","sample")
colors=c("PLF"="#D1832C","pre-PLF"="#44813B","EMT"="#5B93D4","pre-EMT"="#529488","TSL"="#BD3934","NLT"="#A97598","unrecognized"="#C7B299")

c$subtype=factor(c$subtype,levels=c("PLF","pre-PLF","EMT","pre-EMT","TSL","NLT","unrecognized"))

p<-ggplot(c,aes(x=sample,y=percent,fill=subtype))+
	geom_bar(position = "fill",stat="identity")+
	scale_fill_manual(values=colors)

pdf('/public/workspace/liuqzh/gastric/figure2/f2a_A.pdf',width=13,height=4,useDingbats=F)
print(p)
dev.off()
}

#############拟时figure2_b#############
{
bytlib load R-4.0.2
bytlib load gcc
R

library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)

STAD_P <- readRDS("/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/STAD_P_E_T_N.rds")
DefaultAssay(STAD_P)<-'RNA'

STAD_P=subset(STAD_P,cells=which(STAD_P$subtype=='PLF'|STAD_P$subtype=='pre-PLF'|
								 STAD_P$subtype=='EMT'|STAD_P$subtype=='pre-EMT'|
								 STAD_P$subtype=='TSL'|STAD_P$subtype=='NLT'))

STAD_Harmony_Epithelium=STAD_P
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
clustering_DEG_genes <- differentialGeneTest(dat[expressed_genes,], fullModelFormulaStr = '~subtype', cores = 16)
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1401]

#saveRDS(ordering_genes, "/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/mycds_STAD_P_E_T_N_cell_ordering_genes.rds")
#ordering_genes<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/mycds_STAD_P_E_T_N_cell_ordering_genes.rds')

dat <- setOrderingFilter(dat, ordering_genes = ordering_genes)
dat <- reduceDimension(dat, method = 'DDRTree')

dat <- orderCells(dat,root_state = 3)

mycds=dat
saveRDS(mycds, "/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/mycds_STAD_P_E_T_N_cell_SeuratFeatures.rds")
mycds<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/mycds_STAD_P_E_T_N_cell_SeuratFeatures.rds')

colors=c("PLF"="#D1832C","pre-PLF"="#44813B","EMT"="#5B93D4","pre-EMT"="#529488","TSL"="#BD3934","NLT"="#A97598")

##State轨迹分布图
pdf('/public/workspace/liuqzh/gastric_cancer/monocle/GSE150290_monocle_subtype.pdf',width=6,height=6,useDingbats=F)
plot_cell_trajectory(mycds, color_by="subtype",cell_size = 1)+scale_colour_manual(values=colors)
dev.off()

pdf('/public/workspace/liuqzh/gastric_cancer/monocle/GSE150290_monocle_Pseudotime.pdf',width=6,height=6,useDingbats=F)
plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 1)
dev.off()

plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 0.8)
plot_cell_trajectory(mycds, color_by = "Early_Advanced",cell_size = 0.8)

##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters",cell_size = 0.3)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 0.3)
plot4 <- plot_cell_trajectory(mycds, color_by = "celltype",cell_size = 0.3)
plot5 <- plot_cell_trajectory(mycds, color_by = "Early_Advanced",cell_size = 0.1)

pdf('/public/workspace/liuqzh/gastric_cancer/monocle/GSE150290_monocle_state_subtype.pdf',width=18,height=5,useDingbats=F)
plot_cell_trajectory(mycds, color_by = "subtype",cell_size = 1)+facet_wrap(~subtype,nrow = 1)+scale_colour_manual(values=colors)
dev.off()

##合并作图
plotc <- plot1|plot3|plot5
plotc <- plot1|plot5

}

#############HSMM_figure2##############
{
bytlib load R-4.0.2
bytlib load gcc
R

library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

mycds<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/mycds_STAD_P_E_T_N_cell_SeuratFeatures.rds')
ordering_genes<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/mycds_STAD_P_E_T_N_cell_ordering_genes.rds')

colors=c("PLF"="#D1832C","pre-PLF"="#44813B","EMT"="#5B93D4","pre-EMT"="#529488","TSL"="#BD3934","NLT"="#A97598","unrecognized"="#C7B299")

plot_cell_trajectory(mycds, color_by = "subtype",cell_size = 0.8)+scale_colour_manual(values=colors)
#########细胞密度图##########
library(ggpubr)
df=pData(mycds)

p = ggplot(df,aes(Pseudotime,colour=subtype,fill=subtype))+
			geom_density(bw=0.5,size=1,alpha=0.5)+
			theme_classic2()+
			scale_fill_manual(name = "",values = colors)+
			scale_color_manual(name = "",values = colors)

pdf('/public/workspace/liuqzh/gastric_cancer/monocle/GSE150290_monocle_density.pdf',width=6,height=3,useDingbats=F)
print(p)
dev.off()
}

########top5########
bytlib load R-4.0.2
bytlib load gcc
R

library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)
library(patchwork)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(parallel)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(ggpubr)
library(rlang)

mycds<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/mycds_STAD_P_E_T_N_cell_SeuratFeatures.rds')
monocle_meta_data=mycds@phenoData@data
STAD_P <- readRDS("/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/STAD_P_E_T_N.rds")
DefaultAssay(STAD_P)<-'RNA'

STAD_P=subset(STAD_P,cells=which(STAD_P$subtype=='PLF'|STAD_P$subtype=='pre-PLF'|
								 STAD_P$subtype=='EMT'|STAD_P$subtype=='pre-EMT'|
								 STAD_P$subtype=='TSL'|STAD_P$subtype=='NLT'))
STAD_P_RNA <- as.matrix(GetAssayData(STAD_P[["RNA"]], slot='data'))
MP_subtype_III=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv',header=T,sep=',')

#######
MP_1_list=MP_subtype_III$MP_1
STAD_P_RNA_MP1=as.matrix(t(STAD_P_RNA[which(rownames(STAD_P_RNA) %in% MP_1_list),]))

STAD_P_RNA_MP1=cbind(STAD_P_RNA_MP1,monocle_meta_data$Pseudotime)
colnames(STAD_P_RNA_MP1)[101]='Pseudotime'
STAD_P_RNA_MP1=STAD_P_RNA_MP1[which(rownames(STAD_P_RNA_MP1) %in% rownames(STAD_P@meta.data[which(STAD_P@meta.data$subtype %in% c('NLT','TSL','pre-PLF','PLF')),])),]

M_MP1=cor(STAD_P_RNA_MP1,method='pearson')
Pseudotime_MP1=M_MP1[which(colnames(M_MP1) == 'Pseudotime'),-101]
Pseudotime_MP1_list=names(which(Pseudotime_MP1 >= 0.5))

 [1] "ENO1"    "STMN1"   "SNRPG"   "HSPD1"   "H2AFZ"   "HMGB2"   "TUBB"
 [8] "CYCS"    "PLP2"    "SLC25A5" "FABP5"   "CYC1"    "TUBA1B"  "TUBA1C"
[15] "PA2G4"   "SNRPF"   "RAN"     "PSME2"   "NME1"    "SNRPD1"  "SNRPB"
[22] "PSMA7"   "RANBP1"

#######
MP_2_list=MP_subtype_III$MP_2
STAD_P_RNA_MP2=as.matrix(t(STAD_P_RNA[which(rownames(STAD_P_RNA) %in% MP_2_list),]))

STAD_P_RNA_MP2=cbind(STAD_P_RNA_MP2,monocle_meta_data$Pseudotime)
colnames(STAD_P_RNA_MP2)[101]='Pseudotime'
STAD_P_RNA_MP2=STAD_P_RNA_MP2[which(rownames(STAD_P_RNA_MP2) %in% rownames(STAD_P@meta.data[which(STAD_P@meta.data$subtype %in% c('NLT','TSL','pre-EMT','EMT')),])),]

M_MP2=cor(STAD_P_RNA_MP2,method='pearson')
Pseudotime_MP2=M_MP2[which(colnames(M_MP2) == 'Pseudotime'),-101]
Pseudotime_MP2_list=names(which(Pseudotime_MP2 >= 0.5))

 [1] "SFN"      "TMEM54"   "RHOC"     "S100A16"  "S100A14"  "MGST3"
 [7] "MALL"     "SPATS2L"  "MUC13"    "PLAC8"    "CLTB"     "EZR"
[13] "AKR1B10"  "AGPAT2"   "LMO7"     "LGALS3"   "PHGR1"    "C15orf48"
[19] "ISG20"    "IL32"     "CES2"     "CLDN7"    "MYL12B"   "SDCBP2"
[25] "MISP"     "LGALS4"   "CEACAM5"

#######
MP_3_list=MP_subtype_III$MP_3
STAD_P_RNA_MP3=as.matrix(t(STAD_P_RNA[which(rownames(STAD_P_RNA) %in% MP_3_list),]))

STAD_P_RNA_MP3=cbind(STAD_P_RNA_MP3,monocle_meta_data$Pseudotime)
colnames(STAD_P_RNA_MP3)[101]='Pseudotime'

M_MP3=cor(STAD_P_RNA_MP3,method='pearson')
Pseudotime_MP3=M_MP3[which(colnames(M_MP3) == 'Pseudotime'),-101]
Pseudotime_MP3_list=names(which(Pseudotime_MP3 <= -0.5))

[1] "GLUL" "MUC6"

######PLF
colors=c("PLF"="#D1832C","pre-PLF"="#44813B","EMT"="#5B93D4","pre-EMT"="#529488","TSL"="#BD3934","NLT"="#A97598","unrecognized"="#C7B299")

keygenes_MP1=c('H2AFZ','HMGB2','SLC25A5','FABP5','TUBA1B','ENO1')
cds_subset_MP1_all = mycds[keygenes_MP1,which(mycds@phenoData@data$subtype %in% c('NLT','TSL','pre-PLF','PLF'))]

for(i in keygenes_MP1){
	y_names=as.name(i)
	cds_subset_MP1 = cds_subset_MP1_all[i,]
	p1 = plot_genes_in_pseudotime(cds_subset_MP1,color_by = "subtype")+scale_color_manual(name = "",values = colors)

pdf(paste('/public/workspace/liuqzh/gastric_cancer/monocle/',{{y_names}},'_monocle.pdf',sep=''),height=4,width=6)
print(p1)
dev.off()
}

######EMT
colors=c("PLF"="#D1832C","pre-PLF"="#44813B","EMT"="#5B93D4","pre-EMT"="#529488","TSL"="#BD3934","NLT"="#A97598","unrecognized"="#C7B299")

keygenes_MP1=c('S100A16','S100A14','MUC13','LMO7','LGALS3','ISG20','IL32','CLDN7','CEACAM5','C15orf48','TMEM54','MISP')
cds_subset_MP1_all = mycds[keygenes_MP1,which(mycds@phenoData@data$subtype %in% c('NLT','TSL','pre-EMT','EMT'))]

for(i in keygenes_MP1){
	y_names=as.name(i)
	cds_subset_MP1 = cds_subset_MP1_all[i,]
	p1 = plot_genes_in_pseudotime(cds_subset_MP1,color_by = "subtype")+scale_color_manual(name = "",values = colors)

pdf(paste('/public/workspace/liuqzh/gastric_cancer/monocle/',{{y_names}},'_monocle.pdf',sep=''),height=4,width=6)
print(p1)
dev.off()
}

######NLT
colors=c("PLF"="#D1832C","pre-PLF"="#44813B","EMT"="#5B93D4","pre-EMT"="#529488","TSL"="#BD3934","NLT"="#A97598","unrecognized"="#C7B299")

keygenes_MP1=c('GLUL','MUC6')
cds_subset_MP1_all = mycds[keygenes_MP1,which(mycds@phenoData@data$subtype %in% c('NLT','TSL','pre-PLF','PLF','pre-EMT','EMT'))]

for(i in keygenes_MP1){
	y_names=as.name(i)
	cds_subset_MP1 = cds_subset_MP1_all[i,]
	p1 = plot_genes_in_pseudotime(cds_subset_MP1,color_by = "subtype")+scale_color_manual(name = "",values = colors)

pdf(paste('/public/workspace/liuqzh/gastric_cancer/monocle/',{{y_names}},'_monocle.pdf',sep=''),height=4,width=6)
print(p1)
dev.off()
}

