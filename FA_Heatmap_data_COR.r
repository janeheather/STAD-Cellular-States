
########TCGA
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

dat_exp<-read.table('/public/workspace/liuqzh/gastric_cancer/cibersort/STAD_FPKM.txt',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

Hallmarks<-readRDS('/public/workspace/liuqzh/gastric_cancer/geneList_msigdb_Hallmarks.RDS')
Hallmarks_names=c('HALLMARK_FATTY_ACID_METABOLISM')
Hallmarks_gene=Hallmarks[[Hallmarks_names]]

key_gene=c('CD36','FABP1','FABP2','FASN','PPARA','GLUT1','ACC','ACLY','ACSS1','ACSS2','H2AFZ','E2F2')

result_corr=dat_exp[which(rownames(dat_exp) %in% key_gene),]
result_corr=as.data.frame(t(result_corr))

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker=result_Prolif_marker[,'E2F2',drop=F]

result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value[,'E2F2',drop=F]
Heatmap_data=cbind(result_Prolif_marker,data_p_value)
colnames(Heatmap_data)=c('TCGA_cor','TCGA_p_value')

#write.table(Heatmap_data, "/public/workspace/liuqzh/gastric_cancer/bulk_result/TCGA_FA_Heatmap_data_COR.txt",row.names=T,col.names=T,quote=F,sep="\t")
}

########GSE62254
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

key_gene=c('CD36','FABP1','FABP2','FASN','PPARA','GLUT1','ACC','ACLY','ACSS1','ACSS2','H2AFZ','E2F2')

result_corr=dat_exp[which(rownames(dat_exp) %in% key_gene),]
result_corr=as.data.frame(t(result_corr))

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker=result_Prolif_marker[,'E2F2',drop=F]

result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value[,'E2F2',drop=F]
Heatmap_data=cbind(result_Prolif_marker,data_p_value)
colnames(Heatmap_data)=c('GSE62254_cor','GSE62254_p_value')

#write.table(Heatmap_data, "/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE62254_FA_Heatmap_data_COR.txt",row.names=T,col.names=T,quote=F,sep="\t")
}

########GSE26253
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
key_gene=c('CD36','FABP1','FABP2','FASN','PPARA','GLUT1','ACC','ACLY','ACSS1','ACSS2','H2AFZ','E2F2')

result_corr=dat_exp[which(rownames(dat_exp) %in% key_gene),]
result_corr=as.data.frame(t(result_corr))

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker=result_Prolif_marker[,'E2F2',drop=F]

result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value[,'E2F2',drop=F]
Heatmap_data=cbind(result_Prolif_marker,data_p_value)
colnames(Heatmap_data)=c('GSE26253_cor','GSE26253_p_value')

#write.table(Heatmap_data, "/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE26253_FA_Heatmap_data_COR.txt",row.names=T,col.names=T,quote=F,sep="\t")
}

########GSE84426
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

key_gene=c('CD36','FABP1','FABP2','FASN','PPARA','GLUT1','ACC','ACLY','ACSS1','ACSS2','H2AFZ','E2F2')

result_corr=dat_exp[which(rownames(dat_exp) %in% key_gene),]
result_corr=as.data.frame(t(result_corr))

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker=result_Prolif_marker[,'E2F2',drop=F]

result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value[,'E2F2',drop=F]
Heatmap_data=cbind(result_Prolif_marker,data_p_value)
colnames(Heatmap_data)=c('GSE84426_cor','GSE84426_p_value')

#write.table(Heatmap_data, "/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE84426_FA_Heatmap_data_COR.txt",row.names=T,col.names=T,quote=F,sep="\t")
}

########GSE15459
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

key_gene=c('CD36','FABP1','FABP2','FASN','PPARA','GLUT1','ACC','ACLY','ACSS1','ACSS2','H2AFZ','E2F2')

result_corr=dat_exp[which(rownames(dat_exp) %in% key_gene),]
result_corr=as.data.frame(t(result_corr))

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker=result_Prolif_marker[,'E2F2',drop=F]

result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value[,'E2F2',drop=F]
Heatmap_data=cbind(result_Prolif_marker,data_p_value)
colnames(Heatmap_data)=c('GSE15459_cor','GSE15459_p_value')

#write.table(Heatmap_data, "/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE15459_FA_Heatmap_data_COR.txt",row.names=T,col.names=T,quote=F,sep="\t")
}

########GSE113255
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

set.seed(12345)
GSE113255_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE113255/GSE113255_data.rds')
dat_exp=GSE113255_data
dat_exp<-dat_exp[rowSums(dat_exp)>0,]
dat_exp = as.matrix(dat_exp)

key_gene=c('CD36','FABP1','FABP2','FASN','PPARA','GLUT1','ACC','ACLY','ACSS1','ACSS2','H2AFZ','E2F2')

result_corr=dat_exp[which(rownames(dat_exp) %in% key_gene),]
result_corr=as.data.frame(t(result_corr))

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker=result_Prolif_marker[,'E2F2',drop=F]

result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value[,'E2F2',drop=F]
Heatmap_data=cbind(result_Prolif_marker,data_p_value)
colnames(Heatmap_data)=c('GSE113255_cor','GSE113255_p_value')
Heatmap_data=Heatmap_data[sort(rownames(Heatmap_data)),]

#write.table(Heatmap_data, "/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE113255_FA_Heatmap_data_COR.txt",row.names=T,col.names=T,quote=F,sep="\t")
}

########GSE54129
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

set.seed(12345)
GSE54129_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE54129/GSE54129_data.rds')
dat_exp=GSE54129_data
dat_exp<-dat_exp[rowSums(dat_exp)>0,]
dat_exp = as.matrix(dat_exp)

key_gene=c('CD36','FABP1','FABP2','FASN','PPARA','GLUT1','ACC','ACLY','ACSS1','ACSS2','H2AFZ','E2F2')

result_corr=dat_exp[which(rownames(dat_exp) %in% key_gene),]
result_corr=as.data.frame(t(result_corr))

result_Prolif_marker=cor(result_corr,method='pearson')
result_Prolif_marker=result_Prolif_marker[,'E2F2',drop=F]

result_Prolif_marker_p_value=cor.mtest(result_corr,method='pearson')
result_Prolif_marker_p_value=result_Prolif_marker_p_value$p

data_p_value=result_Prolif_marker_p_value[,'E2F2',drop=F]
Heatmap_data=cbind(result_Prolif_marker,data_p_value)
colnames(Heatmap_data)=c('GSE54129_cor','GSE54129_p_value')

#write.table(Heatmap_data, "/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE54129_FA_Heatmap_data_COR.txt",row.names=T,col.names=T,quote=F,sep="\t")
}

#################
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(pheatmap)

###
TCGA_FA_Heatmap_data_COR=read.table('/public/workspace/liuqzh/gastric_cancer/bulk_result/TCGA_FA_Heatmap_data_COR.txt',header=T,sep='\t')
###
GSE84426_FA_Heatmap_data_COR=read.table('/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE84426_FA_Heatmap_data_COR.txt',header=T,sep='\t')
###
GSE62254_FA_Heatmap_data_COR=read.table('/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE62254_FA_Heatmap_data_COR.txt',header=T,sep='\t')
###
GSE54129_FA_Heatmap_data_COR=read.table('/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE54129_FA_Heatmap_data_COR.txt',header=T,sep='\t')
###
GSE26253_FA_Heatmap_data_COR=read.table('/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE26253_FA_Heatmap_data_COR.txt',header=T,sep='\t')
###
GSE15459_FA_Heatmap_data_COR=read.table('/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE15459_FA_Heatmap_data_COR.txt',header=T,sep='\t')
###
GSE113255_FA_Heatmap_data_COR=read.table('/public/workspace/liuqzh/gastric_cancer/bulk_result/GSE113255_FA_Heatmap_data_COR.txt',header=T,sep='\t')
###

heatmap_data= cbind(TCGA_FA_Heatmap_data_COR$TCGA_cor,
					GSE84426_FA_Heatmap_data_COR$GSE84426_cor,
					GSE62254_FA_Heatmap_data_COR$GSE62254_cor,
					GSE54129_FA_Heatmap_data_COR$GSE54129_cor,
					GSE26253_FA_Heatmap_data_COR$GSE26253_cor,
					GSE15459_FA_Heatmap_data_COR$GSE15459_cor,
					GSE113255_FA_Heatmap_data_COR$GSE113255_cor)
	
colnames(heatmap_data)=c('TCGA_cor',
						 'GSE84426_cor',
						 'GSE62254_cor',
						 'GSE54129_cor',
						 'GSE26253_cor',
						 'GSE15459_cor',
						 'GSE113255_cor')

rownames(heatmap_data)=rownames(TCGA_FA_Heatmap_data_COR)
cmt=heatmap_data[-5,]
########
p_value_data= cbind(TCGA_FA_Heatmap_data_COR$TCGA_p_value,
					GSE84426_FA_Heatmap_data_COR$GSE84426_p_value,
					GSE62254_FA_Heatmap_data_COR$GSE62254_p_value,
					GSE54129_FA_Heatmap_data_COR$GSE54129_p_value,
					GSE26253_FA_Heatmap_data_COR$GSE26253_p_value,
					GSE15459_FA_Heatmap_data_COR$GSE15459_p_value,
					GSE113255_FA_Heatmap_data_COR$GSE113255_p_value)
	
colnames(p_value_data)=c('TCGA_p_value',
						 'GSE84426_p_value',
						 'GSE62254_p_value',
						 'GSE54129_p_value',
						 'GSE26253_p_value',
						 'GSE15459_p_value',
						 'GSE113255_p_value')

rownames(p_value_data)=rownames(TCGA_FA_Heatmap_data_COR)

#判断显著性
pmt=p_value_data[-5,]
if (!is.null(pmt)){
	ssmt <- pmt< 0.01
	pmt[ssmt] <-'**'
	smt <- pmt >0.01& pmt <0.05
	pmt[smt] <- '*'
	pmt[!ssmt&!smt]<- ''
	} else {
		pmt <- F
}

#可视化
pdf('/public/workspace/liuqzh/gastric_cancer/bulk_result/FA_Heatmap_data_COR.pdf',height=5,width=4)
pheatmap(cmt,scale = "none",cluster_row = F, cluster_col = F, border=NA,
      display_numbers = pmt,fontsize_number = 12, number_color = "white",
      cellwidth = 20, cellheight =20,
	  color = colorRampPalette(c("#377EB8","white","#E41A1C"))(100))
dev.off()


