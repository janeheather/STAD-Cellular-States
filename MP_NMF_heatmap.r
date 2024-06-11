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

Idents(dat)=dat$Subtype
sig_marker=FindAllMarkers(object = dat, only.pos = TRUE, logfc.threshold = 0.3)

#write.csv(sig_marker, "/public/workspace/liuqzh/gastric_cancer/sig_marker.csv")

#############
bytlib load R-4.0.2
bytlib load gcc
R

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(pheatmap)
DAT<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(DAT) <- "RNA"

table(DAT$Subtype)
a<-table(DAT$orig.ident) %>% as.data.frame()
b<-a[,1]
c<-data.frame()
for(i in b){
	tmp<-subset(DAT,subset=orig.ident==i)
	tmp1<-table(tmp$Subtype) %>% as.data.frame()
	tmp2<-cbind(tmp1,i)
	c<-rbind(c,tmp2)
}
colnames(c)<-c("Subtype","percent","sample")
colors=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")

c$Subtype=factor(c$Subtype,levels=c("MP_1","MP_2","MP_3"))
#write.table(c, "/public/workspace/liuqzh/gastric_cancer/MP1_MP2_MP3_cellmarker.txt", row.names = T, sep = "\t")

c=read.table("/public/workspace/liuqzh/gastric_cancer/MP1_MP2_MP3_cellmarker.txt",header=T,sep="\t")
######
pct = c %>%
	group_by(sample) %>%
	mutate(percent = percent/sum(percent))
	
pct=as.data.frame(pct)

df<-spread(data=pct,key=sample,value=percent)
rownames(df)=df[,1]
df=df[,-c(1,54)]

df=as.data.frame(df)
#write.table(df, "/public/workspace/liuqzh/gastric_cancer/MP1_MP2_MP3_df.txt", row.names = T, sep = "\t")
df=read.table("/public/workspace/liuqzh/gastric_cancer/MP_NMF_heatmap.txt",header=T,sep="\t",row.names=1)
######
library(RColorBrewer)
library(viridis)

pdf('/public/workspace/liuqzh/gastric_cancer/ste_score/MP_NMF_heatmap.pdf',width=12,height=3)
pheatmap(df, cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#5B93D4", "white", "#BD3934"))(256))
dev.off()

###########
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
dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(dat) <- "RNA"

dat$type<-'unrecognized'
dat$type[which(dat$celltype %in% c('B cell','Macrophage','Mast cell','NK cell','T cell'))]<-'Immune'
dat$type[which(dat$celltype %in% c('Endothelial','Fibroblast'))]<-'Stromal'
dat$type[which(dat$celltype %in% c('Epithelium'))]<-'Tumor'
dat<-subset(dat,subset=type %in% c('Immune','Stromal','Tumor'))

DAT=dat
table(DAT$type)
a<-table(DAT$orig.ident) %>% as.data.frame()
b<-a[,1]
c<-data.frame()
for(i in b){
	tmp<-subset(DAT,subset=orig.ident==i)
	tmp1<-table(tmp$type) %>% as.data.frame()
	tmp2<-cbind(tmp1,i)
	c<-rbind(c,tmp2)
}
colnames(c)<-c("Subtype","percent","sample")
colors=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")


