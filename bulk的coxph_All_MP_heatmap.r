
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
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype.csv",sep = ",",row.names = 1,header = F)
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

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)

sur$PATIENT=paste(sur$PATIENT,'-01A',sep='')
ssGSEA_Score$PATIENT=rownames(ssGSEA_Score)

sur_data=merge(ssGSEA_Score,sur,'PATIENT')
rownames(sur_data)=sur_data$PATIENT
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(34,35)])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(34,35)]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(34,35)])

#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_TCGA.rds')
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
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype.csv",sep = ",",row.names = 1,header = F)
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

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=sur$GEO_ID

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

sur_data=sur_data[,-c(34:38)]
sur_data=sur_data[,-c(36:40)]

result=NULL
for(i in colnames(sur_data[,-c(34,35)])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(34,35)]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(34,35)])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE62254.rds')
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
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype.csv",sep = ",",row.names = 1,header = F)
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

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

sur_data=sur_data[,-c(34:41)]
sur_data=sur_data[,-c(36:38)]

result=NULL
for(i in colnames(sur_data[,-c(34,35)])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(34,35)]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(34,35)])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE13861.rds')
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
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype.csv",sep = ",",row.names = 1,header = F)
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

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

sur_data=sur_data[,-c(34:40)]
sur_data=sur_data[,-c(36:38)]

result=NULL
for(i in colnames(sur_data[,-c(34,35)])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(34,35)]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(34,35)])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE26253.rds')
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
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype.csv",sep = ",",row.names = 1,header = F)
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

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur<-read.csv('/public/workspace/yumiao/STAD/bulk/GSE84426_infor.csv',header=T,row.names=1)
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

result=NULL
for(i in colnames(sur_data[,-c(34,35)])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(34,35)]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(34,35)])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE84426.rds')
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
GSE15459_info=read.csv('/public/workspace/liuqzh/STAD_bulk_DATA/GSE15459_data/GSE15459_data_clinical_information.csv',row.names=1)
GSE15459_info$sample=rownames(GSE15459_info)

set.seed(520)
TPM=GSE15459_data
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype.csv",sep = ",",row.names = 1,header = F)
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

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur<-GSE15459_info

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]
sur_data$OS_time=sur_data$OS
sur_data=sur_data[,-c(34:39)]

result=NULL
for(i in colnames(sur_data[,-c(34,35)])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(34,35)]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(34,35)])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE15459.rds')
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
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype.csv",sep = ",",row.names = 1,header = F)
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

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=sur$GEO_ID

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

sur_data=sur_data[,-c(34:36)]
sur_data=sur_data[,-c(36:42)]

result=NULL
for(i in colnames(sur_data[,-c(34,35)])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(34,35)]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(34,35)])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE62254_DFS.rds')
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
genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype.csv",sep = ",",row.names = 1,header = F)
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

scdata = scrabble::score(TPM,
                   groups=geneset,
                   binmat = NULL,
                   bins = NULL,
                   controls = NULL,
                   bin.control = F,
                   center = T,
                   nbin = 30,
                   n = 100,
                   replace = T)

ssGSEA_Score=scdata
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=rownames(sur)

sur_data=merge(ssGSEA_Score,sur,by='sample')
rownames(sur_data)=sur_data$sample
sur_data=sur_data[,-1]

sur_data=sur_data[,-c(34:43)]
sur_data=sur_data[,-c(36)]

result=NULL
for(i in colnames(sur_data[,-c(34,35)])){
	cox_data=as.data.frame(cbind(sur_data[[i]],sur_data[,c(34,35)]))
	names(cox_data)=c('Value','OS','OS.time')
	model <- coxph(Surv(OS.time, OS) ~ Value, data = cox_data)
	coeff=model$coefficients
	p_value=as.numeric(broom::tidy(model)[5])
	data_coxph=cbind(coeff,p_value)
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(34,35)])
#saveRDS(result,'/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE13861_DFS.rds')
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
coxph_TCGA=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_TCGA.rds')
coxph_TCGA=as.data.frame(coxph_TCGA)
colnames(coxph_TCGA)=c('TCGA_coeff','TCGA_p_value')
###
coxph_GSE62254=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE62254.rds')
coxph_GSE62254=as.data.frame(coxph_GSE62254)
colnames(coxph_GSE62254)=c('GSE62254_coeff','GSE62254_p_value')
###
coxph_GSE13861=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE13861.rds')
coxph_GSE13861=as.data.frame(coxph_GSE13861)
colnames(coxph_GSE13861)=c('GSE13861_coeff','GSE13861_p_value')
###
coxph_GSE26253=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE26253.rds')
coxph_GSE26253=as.data.frame(coxph_GSE26253)
colnames(coxph_GSE26253)=c('GSE26253_coeff','GSE26253_p_value')
###
coxph_GSE84426=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE84426.rds')
coxph_GSE84426=as.data.frame(coxph_GSE84426)
colnames(coxph_GSE84426)=c('GSE84426_coeff','GSE84426_p_value')
###
coxph_GSE15459=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE15459.rds')
coxph_GSE15459=as.data.frame(coxph_GSE15459)
colnames(coxph_GSE15459)=c('GSE15459_coeff','GSE15459_p_value')
###
coxph_GSE62254_DFS=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE62254_DFS.rds')
coxph_GSE62254_DFS=as.data.frame(coxph_GSE62254_DFS)
colnames(coxph_GSE62254_DFS)=c('GSE62254_DFS_coeff','GSE62254_DFS_p_value')
###
coxph_GSE13861_DFS=readRDS('/public/workspace/liuqzh/gastric_cancer/figure3/coxph_GSE13861_DFS.rds')
coxph_GSE13861_DFS=as.data.frame(coxph_GSE13861_DFS)
colnames(coxph_GSE13861_DFS)=c('GSE13861_DFS_coeff','GSE13861_DFS_p_value')
###

heatmap_data= cbind(coxph_TCGA$TCGA_coeff,coxph_GSE62254$GSE62254_coeff,
					coxph_GSE13861$GSE13861_coeff,coxph_GSE26253$GSE26253_coeff,
					coxph_GSE84426$GSE84426_coeff,coxph_GSE15459$GSE15459_coeff,
					coxph_GSE62254_DFS$GSE62254_DFS_coeff,coxph_GSE13861_DFS$GSE13861_DFS_coeff)
	
colnames(heatmap_data)=c('TCGA_coeff','GSE62254_coeff',
						 'GSE13861_coeff','GSE26253_coeff',
						 'GSE84426_coeff','GSE15459_coeff',
						 'GSE62254_DFS_coeff','GSE13861_DFS_coeff')

rownames(heatmap_data)=rownames(coxph_TCGA)
cmt=heatmap_data
########
p_value_data= cbind(coxph_TCGA$TCGA_p_value,coxph_GSE62254$GSE62254_p_value,
					coxph_GSE13861$GSE13861_p_value,coxph_GSE26253$GSE26253_p_value,
					coxph_GSE84426$GSE84426_p_value,coxph_GSE15459$GSE15459_p_value,
					coxph_GSE62254_DFS$GSE62254_DFS_p_value,coxph_GSE13861_DFS$GSE13861_DFS_p_value)
	
colnames(p_value_data)=c('TCGA_p_value','GSE62254_p_value',
						 'GSE13861_p_value','GSE26253_p_value',
						 'GSE84426_p_value','GSE15459_p_value',
						 'GSE62254_DFS_p_value','GSE13861_DFS_p_value')

rownames(p_value_data)=rownames(coxph_TCGA)

#判断显著性
pmt=p_value_data
if (!is.null(pmt)){
	ssmt <- pmt< 0.01
	pmt[ssmt] <-'**'
	smt <- pmt >0.01& pmt <0.05
	pmt[smt] <- '*'
	pmt[!ssmt&!smt]<- ''
	} else {
		pmt <- F
}

heatmap_cluster=read.table('/public/workspace/liuqzh/gastric_cancer/figure3/heatmap_cluster.txt',header=T,sep='\t',row.names=1)

cmt=cmt[intersect(rownames(heatmap_cluster),rownames(cmt)),]
pmt=pmt[intersect(rownames(heatmap_cluster),rownames(pmt)),]

#可视化
pdf('/public/workspace/liuqzh/gastric_cancer/figure3/Heatmap_subtype_coxph.pdf',height=15,width=4)
pheatmap(cmt,scale = "none",cluster_row = F, cluster_col = F, border=NA,
      display_numbers = pmt,fontsize_number = 12, number_color = "white",
      cellwidth = 20, cellheight =20,
	  color = colorRampPalette(c("#377EB8","white","#E41A1C"))(100))
dev.off()
