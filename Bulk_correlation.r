####BULK_correlation####
###########################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
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
library(fmsb)

setwd('/public/workspace/liuqzh/gastric/survival/')
dat_exp<-read.table('/public/workspace/liuqzh/gastric_cancer/cibersort/STAD_FPKM.txt',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]
dat_exp <- dat_exp[rowSums(dat_exp)>0,]
dat_exp = as.matrix(dat_exp)

Hallmarks<-readRDS('/public/workspace/liuqzh/gastric_cancer/geneList_msigdb_Hallmarks.RDS')
Hallmarks_names=c('HALLMARK_APICAL_JUNCTION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
				  'HALLMARK_DNA_REPAIR','HALLMARK_E2F_TARGETS',
				  'HALLMARK_APOPTOSIS','HALLMARK_FATTY_ACID_METABOLISM',
				  'HALLMARK_TGF_BETA_SIGNALING','HALLMARK_HYPOXIA',
				  'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_PI3K_AKT_MTOR_SIGNALING',
				  'HALLMARK_GLYCOLYSIS','HALLMARK_TNFA_SIGNALING_VIA_NFKB')

result=NULL
for(i in Hallmarks_names){
	gene <- Hallmarks[i]
	gsva.res<-gsva(dat_exp,gene,method = "ssgsea",parallel.sz = 10)
	tmp_result=as.matrix(t(gsva.res))
	print(i)
	result=cbind(result,tmp_result)
}
Data_MP1<-list(c("ACTG1","AKR1B10","GPX2","H2AFZ","HMGB2","KIAA0101","MUC5AC","STMN1","TFF1","TUBA1B"))
ssGSEA_Score=gsva(dat_exp,Data_MP1, method='ssgsea',parallel.sz = 10)
ssGSEA_Score=as.matrix(t(ssGSEA_Score))#ssGSEA计算
result_corr=cbind(result,ssGSEA_Score)
colnames(result_corr)[ncol(result_corr)]='Prolif_marker'

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value
data_p_value=data_p_value[-c(nrow(data_p_value)),ncol(data_p_value),drop=F]
data_p_value=-log10(data_p_value)
data_p_value=as.data.frame(t(data_p_value))
data <- rbind(rep(30,12) , rep(0,12) , data_p_value)

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_TCGA_STAD_HALLMARKER_PLF.pdf',height=10,width=10)
radarchart(data,axistype=1, 
	#custom polygon
	pcol='#5B93D4',plwd=2, 
	#custom the grid
	cglcol="grey",cglty=1,axislabcol="grey",caxislabels=seq(0,30,5),cglwd=0.8,
	#custom labels
	vlcex=0.8)
dev.off()
}

###########################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
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
library(fmsb)

setwd('/public/workspace/yumiao/STAD/bulk/GSE62254/')
sur<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE62254/GSE62254.csv',header=T)
gse<-read.table('GSE62254_series_matrix (1).txt',header=T,sep='\t')
ref<-read.table('GPL570-55999 (1).txt',header=T,sep='\t')

colnames(ref)[1]<-'ID_REF'
GSE62254<-merge(ref,gse,by='ID_REF')
final<-aggregate(GSE62254[,3:ncol(GSE62254)],list(GSE62254[,2]),mean)
final<-final[-1,]
rownames(final)<-final[,1]
final<-final[,-1]
#final<-log2(final+1)
sur<-arrange(sur,GEO_ID)
dat_exp=final

set.seed(12345)

dat_exp = as.matrix(dat_exp)

Hallmarks<-readRDS('/public/workspace/liuqzh/gastric_cancer/geneList_msigdb_Hallmarks.RDS')
Hallmarks_names=c('HALLMARK_APICAL_JUNCTION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
				  'HALLMARK_DNA_REPAIR','HALLMARK_E2F_TARGETS',
				  'HALLMARK_APOPTOSIS','HALLMARK_FATTY_ACID_METABOLISM',
				  'HALLMARK_TGF_BETA_SIGNALING','HALLMARK_HYPOXIA',
				  'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_PI3K_AKT_MTOR_SIGNALING',
				  'HALLMARK_GLYCOLYSIS','HALLMARK_TNFA_SIGNALING_VIA_NFKB')

result=NULL
for(i in Hallmarks_names){
	gene <- Hallmarks[i]
	gsva.res<-gsva(dat_exp,gene,method = "ssgsea",parallel.sz = 10)
	tmp_result=as.matrix(t(gsva.res))
	print(i)
	result=cbind(result,tmp_result)
}
Data_MP1<-list(c("ACTG1","AKR1B10","GPX2","H2AFZ","HMGB2","KIAA0101","MUC5AC","STMN1","TFF1","TUBA1B"))
ssGSEA_Score=gsva(dat_exp,Data_MP1, method='ssgsea',parallel.sz = 10)
ssGSEA_Score=as.matrix(t(ssGSEA_Score))#ssGSEA计算
result_corr=cbind(result,ssGSEA_Score)
colnames(result_corr)[ncol(result_corr)]='Prolif_marker'

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value
data_p_value=data_p_value[-c(nrow(data_p_value)),ncol(data_p_value),drop=F]
data_p_value=-log10(data_p_value)
data_p_value=as.data.frame(t(data_p_value))
data <- rbind(rep(30,12) , rep(0,12) , data_p_value)

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE62254_HALLMARKER_PLF.pdf',height=10,width=10)
radarchart(data,axistype=1, 
	#custom polygon
	pcol='#5B93D4',plwd=2, 
	#custom the grid
	cglcol="grey",cglty=1,axislabcol="grey",caxislabels=seq(0,30,5),cglwd=0.8,
	#custom labels
	vlcex=0.8)
dev.off()
}

###########################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
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
library(fmsb)

setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE26253/')
sur<-read.csv('GSE26253_SMC_sur.csv',header=T,row.names=2)
gse<-read.table('GSE26253_matrix.txt',header=T,sep='\t')
ref<-read.table('GPL8432-11703.txt',header=T,sep='\t')

colnames(ref)[1]<-'ID_REF'
GSE26253<-merge(ref,gse,by='ID_REF')
final<-aggregate(GSE26253[,3:ncol(GSE26253)],list(GSE26253$Symbol),mean)
rownames(final)<-final[,1]
final<-final[,-1]
#final<-log2(final+1)
final=as.matrix(final)
sur<-arrange(sur,Patients_ID)

set.seed(12345)
dat_exp=final
dat_exp <- dat_exp[rowSums(dat_exp)>0,]
dat_exp = as.matrix(dat_exp)

Hallmarks<-readRDS('/public/workspace/liuqzh/gastric_cancer/geneList_msigdb_Hallmarks.RDS')
Hallmarks_names=c('HALLMARK_APICAL_JUNCTION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
				  'HALLMARK_DNA_REPAIR','HALLMARK_E2F_TARGETS',
				  'HALLMARK_APOPTOSIS','HALLMARK_FATTY_ACID_METABOLISM',
				  'HALLMARK_TGF_BETA_SIGNALING','HALLMARK_HYPOXIA',
				  'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_PI3K_AKT_MTOR_SIGNALING',
				  'HALLMARK_GLYCOLYSIS','HALLMARK_TNFA_SIGNALING_VIA_NFKB')

result=NULL
for(i in Hallmarks_names){
	gene <- Hallmarks[i]
	gsva.res<-gsva(dat_exp,gene,method = "ssgsea",parallel.sz = 10)
	tmp_result=as.matrix(t(gsva.res))
	print(i)
	result=cbind(result,tmp_result)
}
Data_MP1<-list(c("ACTG1","AKR1B10","GPX2","H2AFZ","HMGB2","KIAA0101","MUC5AC","STMN1","TFF1","TUBA1B"))
ssGSEA_Score=gsva(dat_exp,Data_MP1, method='ssgsea',parallel.sz = 10)
ssGSEA_Score=as.matrix(t(ssGSEA_Score))#ssGSEA计算
result_corr=cbind(result,ssGSEA_Score)
colnames(result_corr)[ncol(result_corr)]='Prolif_marker'

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value
data_p_value=data_p_value[-c(nrow(data_p_value)),ncol(data_p_value),drop=F]
data_p_value=-log10(data_p_value)
data_p_value=as.data.frame(t(data_p_value))
data <- rbind(rep(30,12) , rep(0,12) , data_p_value)

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE26253_HALLMARKER_PLF.pdf',height=10,width=10)
radarchart(data,axistype=1, 
	#custom polygon
	pcol='#5B93D4',plwd=2, 
	#custom the grid
	cglcol="grey",cglty=1,axislabcol="grey",caxislabels=seq(0,30,5),cglwd=0.8,
	#custom labels
	vlcex=0.8)
dev.off()
}

###########################
{
bytlib load gcc
bytlib load R-4.0.2
R

library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
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
library(fmsb)
 
cli<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE84426_infor.csv',header=T)
expr<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE84426_matrix.csv',header=T)
anno<-read.table('/public/workspace/yumiao/STAD/bulk/GPL6947-13512.txt',sep='\t',header=T)
colnames(anno)[1]<-'ID_REF'
dat<-merge(anno,expr,by='ID_REF')
final<-aggregate(dat[,3:ncol(dat)],list(dat$Symbol),mean)
final<-final[-1,]
rownames(final)<-final[,1]
final<-final[,-1]

final<-log2(final+1)

set.seed(12345)
dat_exp=final

dat_exp <- dat_exp[rowSums(dat_exp)>0,]
dat_exp = as.matrix(dat_exp)

Hallmarks<-readRDS('/public/workspace/liuqzh/gastric_cancer/geneList_msigdb_Hallmarks.RDS')
Hallmarks_names=c('HALLMARK_APICAL_JUNCTION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
				  'HALLMARK_DNA_REPAIR','HALLMARK_E2F_TARGETS',
				  'HALLMARK_APOPTOSIS','HALLMARK_FATTY_ACID_METABOLISM',
				  'HALLMARK_TGF_BETA_SIGNALING','HALLMARK_HYPOXIA',
				  'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_PI3K_AKT_MTOR_SIGNALING',
				  'HALLMARK_GLYCOLYSIS','HALLMARK_TNFA_SIGNALING_VIA_NFKB')

result=NULL
for(i in Hallmarks_names){
	gene <- Hallmarks[i]
	gsva.res<-gsva(dat_exp,gene,method = "ssgsea",parallel.sz = 10)
	tmp_result=as.matrix(t(gsva.res))
	print(i)
	result=cbind(result,tmp_result)
}
Data_MP1<-list(c("ACTG1","AKR1B10","GPX2","H2AFZ","HMGB2","KIAA0101","MUC5AC","STMN1","TFF1","TUBA1B"))
ssGSEA_Score=gsva(dat_exp,Data_MP1, method='ssgsea',parallel.sz = 10)
ssGSEA_Score=as.matrix(t(ssGSEA_Score))#ssGSEA计算
result_corr=cbind(result,ssGSEA_Score)
colnames(result_corr)[ncol(result_corr)]='Prolif_marker'

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value
data_p_value=data_p_value[-c(nrow(data_p_value)),ncol(data_p_value),drop=F]
data_p_value=-log10(data_p_value)
data_p_value=as.data.frame(t(data_p_value))
data <- rbind(rep(5,12) , rep(0,12) , data_p_value)

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE84426_HALLMARKER_PLF.pdf',height=10,width=10)
radarchart(data,axistype=1, 
	#custom polygon
	pcol='#5B93D4',plwd=2, 
	#custom the grid
	cglcol="grey",cglty=1,axislabcol="grey",caxislabels=seq(0,5,1),cglwd=0.8,
	#custom labels
	vlcex=0.8)
dev.off()
}

#############
{
bytlib load R-4.0.2
bytlib load gcc
R

library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
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
library(fmsb)

GSE15459_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE15459/GSE15459_data.rds')
GSE15459_info=read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE15459/GSE15459_data_clinical_information.csv',row.names=1)
GSE15459_info$sample=rownames(GSE15459_info)
cli=GSE15459_info

set.seed(12345)
dat_exp=GSE15459_data
dat_exp<-dat_exp[rowSums(dat_exp)>0,]
dat_exp = as.matrix(dat_exp)

Hallmarks<-readRDS('/public/workspace/liuqzh/gastric_cancer/geneList_msigdb_Hallmarks.RDS')
Hallmarks_names=c('HALLMARK_APICAL_JUNCTION','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
				  'HALLMARK_DNA_REPAIR','HALLMARK_E2F_TARGETS',
				  'HALLMARK_APOPTOSIS','HALLMARK_FATTY_ACID_METABOLISM',
				  'HALLMARK_TGF_BETA_SIGNALING','HALLMARK_HYPOXIA',
				  'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_PI3K_AKT_MTOR_SIGNALING',
				  'HALLMARK_GLYCOLYSIS','HALLMARK_TNFA_SIGNALING_VIA_NFKB')

result=NULL
for(i in Hallmarks_names){
	gene <- Hallmarks[i]
	gsva.res<-gsva(dat_exp,gene,method = "ssgsea",parallel.sz = 10)
	tmp_result=as.matrix(t(gsva.res))
	print(i)
	result=cbind(result,tmp_result)
}
Data_MP1<-list(c("ACTG1","AKR1B10","GPX2","H2AFZ","HMGB2","KIAA0101","MUC5AC","STMN1","TFF1","TUBA1B"))
ssGSEA_Score=gsva(dat_exp,Data_MP1, method='ssgsea',parallel.sz = 10)
ssGSEA_Score=as.matrix(t(ssGSEA_Score))#ssGSEA计算
result_corr=cbind(result,ssGSEA_Score)
colnames(result_corr)[ncol(result_corr)]='Prolif_marker'

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value
data_p_value=data_p_value[-c(nrow(data_p_value)),ncol(data_p_value),drop=F]
data_p_value=-log10(data_p_value)
data_p_value=as.data.frame(t(data_p_value))
data <- rbind(rep(30,12) , rep(0,12) , data_p_value)

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/cor_GSE15459_HALLMARKER_PLF.pdf',height=10,width=10)
radarchart(data,axistype=1, 
	#custom polygon
	pcol='#5B93D4',plwd=2, 
	#custom the grid
	cglcol="grey",cglty=1,axislabcol="grey",caxislabels=seq(0,30,5),cglwd=0.8,
	#custom labels
	vlcex=0.8)
dev.off()
}

