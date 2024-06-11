
####OS####
###################
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

setwd('/public/workspace/liuqzh/gastric/survival/')
dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
#write.table(dat_exp,file="/public/workspace/liuqzh/gastric_cancer/bulk_data/TCGA/tcga_dat_expr.txt",row.names=T,col.names=T,quote=F,sep="\t")
sur<-read.table('/public/workspace/yumiao/STAD/tcga/TCGA-STAD.survival.tsv',header=T,sep='\t',row.names=1)

colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]
sur=sur[grep("-01A",rownames(sur)),]

sur=sur[intersect(rownames(sur),colnames(dat_exp)),]
dat_exp=dat_exp[,intersect(rownames(sur),colnames(dat_exp))]
dat_exp<-dat_exp[rowSums(dat_exp)>0,]

sur$OS.time=sur$OS.time/30
sur$OS[which(sur$OS.time >= 50)]=0

setwd("/public/workspace/liuqzh/gastric_cancer/bulk_data/TCGA/")
##将准备好的表达谱保存为txt格式，这里是用ncbiid，如果是用genesymbol,改成id="GeneSymbol"即可
#filterCommonGenes(input.f="tcga_dat_expr.txt", output.f="tcga.gct", id="GeneSymbol")
#estimateScore(input.ds="tcga.gct", output.ds="tcga_estimate_score.gct", platform="affymetrix")
#estimate_score <- read.table("tcga_estimate_score.gct", skip = 2, header = TRUE)
#write.csv(estimate_score,"tcga_est.csv",row.names = FALSE)

set.seed(520)
TPM=dat_exp

####
genelist=rownames(TPM)
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)

sur$PATIENT=paste(sur$PATIENT,'-01A',sep='')
ssGSEA_Score$PATIENT=rownames(ssGSEA_Score)

sur_data=merge(ssGSEA_Score,sur,'PATIENT')
rownames(sur_data)=sur_data$PATIENT
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(ncol(sur_data),c(ncol(sur_data)-1))]))
	names(cox_data)=c('Value','OS.time','OS')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])

#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_TCGA_all_genen.rds')
}

###################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

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

set.seed(520)
TPM=dat_exp

####
genelist=rownames(TPM)
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)
ssGSEA_Score$sample=rownames(ssGSEA_Score)

sur=sur[,c(1,6,7)]
colnames(sur)[1]='sample'

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(ncol(sur_data),c(ncol(sur_data)-1))]))
	names(cox_data)=c('Value','OS.time','OS')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE62254_all_gene.rds')
}

###################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/')
sur<-read.csv('GSE13861_YUHS_sur.csv',header=T,row.names=2)
gse<-read.table('GSE13861_matrix.txt',header=T,sep='\t',row.names=1)
ref<-read.table('GPL6884.txt',header=T,sep='\t',row.names=1)

GSE13861=gse[intersect(rownames(ref),rownames(gse)),]

ref$ID_REF=rownames(ref)
GSE13861$ID_REF=rownames(GSE13861)
GSE13861=merge(GSE13861,ref,'ID_REF')

final<-aggregate(GSE13861[,2:c(ncol(GSE13861)-1)],list(GSE13861$Symbol),mean)
rownames(final)<-final[,1]
final<-final[,-1]
#final<-log2(final+1)
final=as.matrix(final)
sur<-arrange(sur,Patients_ID)
dat_exp=final

set.seed(520)
TPM=dat_exp

####
genelist=rownames(TPM)
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)
ssGSEA_Score$sample=rownames(ssGSEA_Score)

sur=sur[,c(10,9)]
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(ncol(sur_data),c(ncol(sur_data)-1))]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE13861_all_gene.rds')
}

###################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

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

dat_exp=final

set.seed(520)
TPM=dat_exp

####
genelist=rownames(TPM)
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)
ssGSEA_Score$sample=rownames(ssGSEA_Score)

sur=sur[,c(9,8)]
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(ncol(sur_data),c(ncol(sur_data)-1))]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE26253_all_gene.rds')
}

###################
{
bytlib load gcc
bytlib load R-4.0.2
R

library(GSVA)
library(survival)
library(survminer)
library(dplyr)
 
sur<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE84426_infor.csv',header=T,row.names=1)
expr<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE84426_matrix.csv',header=T)
anno<-read.table('/public/workspace/yumiao/STAD/bulk/GPL6947-13512.txt',sep='\t',header=T)
colnames(anno)[1]<-'ID_REF'
dat<-merge(anno,expr,by='ID_REF')
final<-aggregate(dat[,3:ncol(dat)],list(dat$Symbol),mean)
final<-final[-1,]
rownames(final)<-final[,1]
final<-final[,-1]

final<-log2(final+1)

dat_exp=final

set.seed(520)
TPM=dat_exp

####
genelist=rownames(TPM)
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)
ssGSEA_Score$sample=rownames(ssGSEA_Score)

sur=sur[,c(2,1)]
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(ncol(sur_data),c(ncol(sur_data)-1))]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE84426_all_gene.rds')
}

###################
{
bytlib load R-4.0.2
bytlib load gcc
R

library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

GSE15459_data=readRDS('/public/workspace/liuqzh/STAD_bulk_DATA/GSE15459_data/GSE15459_data.rds')
sur=read.csv('/public/workspace/liuqzh/STAD_bulk_DATA/GSE15459_data/GSE15459_data_clinical_information.csv',row.names=1)
sur$sample=rownames(sur)

set.seed(520)
TPM=GSE15459_data

####
genelist=rownames(TPM)
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)
ssGSEA_Score$sample=rownames(ssGSEA_Score)

sur=sur[,c(6,7)]
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(ncol(sur_data),c(ncol(sur_data)-1))]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE15459_all_gene.rds')
}


####DFS####
###################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

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

set.seed(520)
TPM=dat_exp

####
genelist=rownames(TPM)
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)
ssGSEA_Score$sample=rownames(ssGSEA_Score)

sur=sur[,c(1,5,4)]
colnames(sur)[1]='sample'

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(ncol(sur_data),c(ncol(sur_data)-1))]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE62254_DFS_all_gene.rds')
}

###################
{
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/')
sur<-read.csv('GSE13861_YUHS_sur.csv',header=T,row.names=2)
gse<-read.table('GSE13861_matrix.txt',header=T,sep='\t',row.names=1)
ref<-read.table('GPL6884.txt',header=T,sep='\t',row.names=1)

GSE13861=gse[intersect(rownames(ref),rownames(gse)),]

ref$ID_REF=rownames(ref)
GSE13861$ID_REF=rownames(GSE13861)
GSE13861=merge(GSE13861,ref,'ID_REF')

final<-aggregate(GSE13861[,2:c(ncol(GSE13861)-1)],list(GSE13861$Symbol),mean)
rownames(final)<-final[,1]
final<-final[,-1]
#final<-log2(final+1)
final=as.matrix(final)
sur<-arrange(sur,Patients_ID)
dat_exp=final

set.seed(520)
TPM=dat_exp

####
genelist=rownames(TPM)
####
ssGSEA_Score=t(TPM[which(rownames(TPM) %in% genelist),])
ssGSEA_Score=as.data.frame(ssGSEA_Score)
ssGSEA_Score$sample=rownames(ssGSEA_Score)

sur=sur[,c(12,11)]
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(ncol(sur_data),c(ncol(sur_data)-1))]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE13861_DFS_all_gene.rds')
}

###################
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
coxph_TCGA=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_TCGA_all_genen.rds')
coxph_TCGA=as.data.frame(coxph_TCGA)
coxph_TCGA=coxph_TCGA[which(coxph_TCGA$p_value < 0.05 & coxph_TCGA$coeff > 0),]

###
coxph_GSE62254=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE62254_all_gene.rds')
coxph_GSE62254=as.data.frame(coxph_GSE62254)
coxph_GSE62254=coxph_GSE62254[which(coxph_GSE62254$p_value < 0.05 & coxph_GSE62254$coeff > 0),]

###
coxph_GSE13861=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE13861_all_gene.rds')
coxph_GSE13861=as.data.frame(coxph_GSE13861)
coxph_GSE13861=coxph_GSE13861[which(coxph_GSE13861$p_value < 0.05 & coxph_GSE13861$coeff > 0),]

###
coxph_GSE26253=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE26253_all_gene.rds')
coxph_GSE26253=as.data.frame(coxph_GSE26253)
coxph_GSE26253=coxph_GSE26253[which(coxph_GSE26253$p_value < 0.05 & coxph_GSE26253$coeff > 0),]

###
coxph_GSE84426=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE84426_all_gene.rds')
coxph_GSE84426=as.data.frame(coxph_GSE84426)
coxph_GSE84426=coxph_GSE84426[which(coxph_GSE84426$p_value < 0.05 & coxph_GSE84426$coeff > 0),]

###
coxph_GSE15459=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE15459_all_gene.rds')
coxph_GSE15459=as.data.frame(coxph_GSE15459)
coxph_GSE15459=coxph_GSE15459[which(coxph_GSE15459$p_value < 0.05 & coxph_GSE15459$coeff > 0),]

###
coxph_GSE62254_DFS=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE62254_DFS_all_gene.rds')
coxph_GSE62254_DFS=as.data.frame(coxph_GSE62254_DFS)
coxph_GSE62254_DFS=coxph_GSE62254_DFS[which(coxph_GSE62254_DFS$p_value < 0.05 & coxph_GSE62254_DFS$coeff > 0),]

###
coxph_GSE13861_DFS=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE13861_DFS_all_gene.rds')
coxph_GSE13861_DFS=as.data.frame(coxph_GSE13861_DFS)
coxph_GSE13861_DFS=coxph_GSE13861_DFS[which(coxph_GSE13861_DFS$p_value < 0.05 & coxph_GSE13861_DFS$coeff > 0),]

###
gene=Reduce(intersect,list(rownames(coxph_TCGA),
					  rownames(coxph_GSE62254),
					  rownames(coxph_GSE13861),
					  rownames(coxph_GSE26253),
					  rownames(coxph_GSE84426),
					  rownames(coxph_GSE15459),
					  rownames(coxph_GSE62254_DFS),
					  rownames(coxph_GSE13861_DFS)))

############
bytlib load R-4.0.2
bytlib load gcc
R
library(Seurat)
library(monocle)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)

set.seed(520)
STAD_Tr_info<-readRDS("/public/workspace/liuqzh/gastric/STAD-PRJEB25780_info.Rds")
#STAD_Tr_info=as.data.frame(STAD_Tr_info[1:45,])
#rownames(STAD_Tr_info)=STAD_Tr_info[,1]

STAD_Tr_data<-readRDS("/public/workspace/liuqzh/gastric/STAD-PRJEB25780.Response.Rds")
STAD_Tr_data=as.data.frame(STAD_Tr_data)
rownames(STAD_Tr_data)=STAD_Tr_data[,1]
STAD_Tr_data=STAD_Tr_data[,-1]


STAD_Tr_data<-t(STAD_Tr_data) %>% as.data.frame
STAD_Tr_data$sample_id<-rownames(STAD_Tr_data)
dat<-merge(STAD_Tr_data,STAD_Tr_info,by='sample_id')
rownames(dat)<-paste0(dat$response_NR,'_',dat$sample_id)
dat<-dat[,-1]
dat<-dat[,1:(ncol(STAD_Tr_data)-1)]
dat<-t(dat) %>% as.data.frame

pvalue<-apply(dat,1,function(x) wilcox.test(x[grep("^R",colnames(dat))],x[grep("^N",colnames(dat))],paired = F)$p.value)
log2FC<-apply(dat,1,function(x){
  log2( (mean(x[grep("^R",colnames(dat))])+1) / (mean(x[grep("^N",colnames(dat))])+1) )
  
  
})

dat$p<-pvalue
dat$log2FC<-log2FC


dat1<-filter(dat,pvalue<0.05,abs(log2FC)>0.8)
rownames(dat1[which(dat1$log2FC < -1),])
