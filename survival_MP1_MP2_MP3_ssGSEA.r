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
sur$OS[which(sur$OS.time >= 60)]=0

setwd("/public/workspace/liuqzh/gastric_cancer/bulk_data/TCGA/")
##将准备好的表达谱保存为txt格式，这里是用ncbiid，如果是用genesymbol,改成id="GeneSymbol"即可
#filterCommonGenes(input.f="tcga_dat_expr.txt", output.f="tcga.gct", id="GeneSymbol")
#estimateScore(input.ds="tcga.gct", output.ds="tcga_estimate_score.gct", platform="affymetrix")
#estimate_score <- read.table("tcga_estimate_score.gct", skip = 2, header = TRUE)
#write.csv(estimate_score,"tcga_est.csv",row.names = FALSE)

set.seed(520)
Data=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv',header=T,sep=',')
Data_MP1<-list(Data$MP_1)
Data_MP2<-list(Data$MP_2)
Data_MP3<-list(Data$MP_3)

ssGSEA_Score_MP1=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP1, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP2=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP2, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP3=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP3, method='ssgsea', kcdf='Poisson')))#ssGSEA计算

colnames(ssGSEA_Score_MP1)='MP1'
colnames(ssGSEA_Score_MP2)='MP2'
colnames(ssGSEA_Score_MP3)='MP3'

ssGSEA_Score=cbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$MP1_Sub='Other'
ssGSEA_Score$MP1_Sub[which(ssGSEA_Score$MP1 > quantile(ssGSEA_Score$MP1,0.5))]='MP1'
ssGSEA_Score$MP2_Sub='Other'
ssGSEA_Score$MP2_Sub[which(ssGSEA_Score$MP2 > quantile(ssGSEA_Score$MP2,0.5))]='MP2'
ssGSEA_Score$MP3_Sub='Other'
ssGSEA_Score$MP3_Sub[which(ssGSEA_Score$MP3 > quantile(ssGSEA_Score$MP3,0.5))]='MP3'

ssGSEA_Score$Subtype_State='Variable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='MP1' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='MP2' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='MP3')]='Stable State'

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=rownames(sur)

cli=merge(sur,ssGSEA_Score,by='sample')
rownames(cli)=cli$sample

cli$group<-factor(cli$Subtype_State,levels=c('Variable State','Stable State'))

fit <- survfit(Surv(OS.time,OS) ~ group, data = cli)

pdf('/public/workspace/liuqzh/gastric_cancer/survival/TCGA_STAD_State.pdf',height=5,width=5)
ggsurvplot(
	fit,
	data = cli,
	title = "Kaplan-Meier Survival Curve",
	pval = TRUE,  
	pval.method = TRUE,  
	pval.size = 6,  
	pval.coord = c(60, 0.1),  
	legend.title = "Group",
	xlab = "Time (Months)",
	ylab = "Survival Probability",
	font.main = 14,  
	font.x = 12,  
	font.y = 12,  
	font.tickslab = 10,  
	palette = c("#377EB8","#E41A1C"),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('TCGA_STAD')
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

Data=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv',header=T,sep=',')
Data_MP1<-list(Data$MP_1)
Data_MP2<-list(Data$MP_2)
Data_MP3<-list(Data$MP_3)

ssGSEA_Score_MP1=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP1, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP2=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP2, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP3=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP3, method='ssgsea', kcdf='Poisson')))#ssGSEA计算

colnames(ssGSEA_Score_MP1)='MP1'
colnames(ssGSEA_Score_MP2)='MP2'
colnames(ssGSEA_Score_MP3)='MP3'

ssGSEA_Score=cbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$MP1_Sub='Other'
ssGSEA_Score$MP1_Sub[which(ssGSEA_Score$MP1 > quantile(ssGSEA_Score$MP1,0.5))]='MP1'
ssGSEA_Score$MP2_Sub='Other'
ssGSEA_Score$MP2_Sub[which(ssGSEA_Score$MP2 > quantile(ssGSEA_Score$MP2,0.5))]='MP2'
ssGSEA_Score$MP3_Sub='Other'
ssGSEA_Score$MP3_Sub[which(ssGSEA_Score$MP3 > quantile(ssGSEA_Score$MP3,0.5))]='MP3'

ssGSEA_Score$Subtype_State='Variable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='MP1' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='MP2' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='MP3')]='Stable State'

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=sur$GEO_ID

cli=merge(sur,ssGSEA_Score,by='sample')
rownames(cli)=cli$sample

cli$group<-factor(cli$Subtype_State,levels=c('Variable State','Stable State'))

fit <- survfit(Surv(OS_m,Death) ~ group, data = cli)

pdf('/public/workspace/liuqzh/gastric_cancer/survival/GSE62254_State.pdf',height=5,width=5)
ggsurvplot(
	fit,
	data = cli,
	title = "Kaplan-Meier Survival Curve",
	pval = TRUE,  
	pval.method = TRUE,  
	pval.size = 6,  
	pval.coord = c(60, 0.1),  
	legend.title = "Group",
	xlab = "Time (Months)",
	ylab = "Survival Probability",
	font.main = 14,  
	font.x = 12,  
	font.y = 12,  
	font.tickslab = 10,  
	palette = c("#377EB8","#E41A1C"),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('GSE62254')
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

setwd('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/')
sur<-read.csv('GSE13861_YUHS_sur.csv',header=T,row.names=2)
gse<-read.table('GSE13861_matrix.txt',header=T,sep='\t',row.names=1)
ref<-read.table('GPL6884.txt',header=T,sep='\t',row.names=1)

GSE13861=gse[intersect(rownames(ref),rownames(gse)),]

ref$ID_REF=rownames(ref)
GSE13861$ID_REF=rownames(GSE13861)
GSE13861=merge(GSE13861,ref,'ID_REF')

final<-aggregate(GSE13861[,2:(ncol(GSE13861)-1)],list(GSE13861$Symbol),mean)
rownames(final)<-final[,1]
final<-final[,-1]
#final<-log2(final+1)
final=as.matrix(final)
sur<-arrange(sur,Patients_ID)
dat_exp=final

Data=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv',header=T,sep=',')
Data_MP1<-list(Data$MP_1)
Data_MP2<-list(Data$MP_2)
Data_MP3<-list(Data$MP_3)

ssGSEA_Score_MP1=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP1, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP2=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP2, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP3=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP3, method='ssgsea', kcdf='Poisson')))#ssGSEA计算

colnames(ssGSEA_Score_MP1)='MP1'
colnames(ssGSEA_Score_MP2)='MP2'
colnames(ssGSEA_Score_MP3)='MP3'

ssGSEA_Score=cbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$MP1_Sub='Other'
ssGSEA_Score$MP1_Sub[which(ssGSEA_Score$MP1 > quantile(ssGSEA_Score$MP1,0.5))]='MP1'
ssGSEA_Score$MP2_Sub='Other'
ssGSEA_Score$MP2_Sub[which(ssGSEA_Score$MP2 > quantile(ssGSEA_Score$MP2,0.5))]='MP2'
ssGSEA_Score$MP3_Sub='Other'
ssGSEA_Score$MP3_Sub[which(ssGSEA_Score$MP3 > quantile(ssGSEA_Score$MP3,0.5))]='MP3'

ssGSEA_Score$Subtype_State='Variable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='MP1' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='MP2' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='MP3')]='Stable State'

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=rownames(sur)

cli=merge(sur,ssGSEA_Score,by='sample')
rownames(cli)=cli$sample

cli$group<-factor(cli$Subtype_State,levels=c('Variable State','Stable State'))

fit <- survfit(Surv(OS_m,Death) ~ group, data = cli)

pdf('/public/workspace/liuqzh/gastric_cancer/survival/GSE13861_State.pdf',height=5,width=5)
ggsurvplot(
	fit,
	data = cli,
	title = "Kaplan-Meier Survival Curve",
	pval = TRUE,  
	pval.method = TRUE,  
	pval.size = 6,  
	pval.coord = c(60, 0.1),  
	legend.title = "Group",
	xlab = "Time (Months)",
	ylab = "Survival Probability",
	font.main = 14,  
	font.x = 12,  
	font.y = 12,  
	font.tickslab = 10,  
	palette = c("#377EB8","#E41A1C"),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('GSE13861')
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

Data=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv',header=T,sep=',')
Data_MP1<-list(Data$MP_1)
Data_MP2<-list(Data$MP_2)
Data_MP3<-list(Data$MP_3)

ssGSEA_Score_MP1=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP1, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP2=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP2, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP3=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP3, method='ssgsea', kcdf='Poisson')))#ssGSEA计算

colnames(ssGSEA_Score_MP1)='MP1'
colnames(ssGSEA_Score_MP2)='MP2'
colnames(ssGSEA_Score_MP3)='MP3'

ssGSEA_Score=cbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$MP1_Sub='Other'
ssGSEA_Score$MP1_Sub[which(ssGSEA_Score$MP1 > quantile(ssGSEA_Score$MP1,0.5))]='MP1'
ssGSEA_Score$MP2_Sub='Other'
ssGSEA_Score$MP2_Sub[which(ssGSEA_Score$MP2 > quantile(ssGSEA_Score$MP2,0.5))]='MP2'
ssGSEA_Score$MP3_Sub='Other'
ssGSEA_Score$MP3_Sub[which(ssGSEA_Score$MP3 > quantile(ssGSEA_Score$MP3,0.5))]='MP3'

ssGSEA_Score$Subtype_State='Variable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='MP1' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='MP2' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='MP3')]='Stable State'

ssGSEA_Score$sample=rownames(ssGSEA_Score)
sur$sample=rownames(sur)

cli=merge(sur,ssGSEA_Score,by='sample')
rownames(cli)=cli$sample

cli$group<-factor(cli$Subtype_State,levels=c('Variable State','Stable State'))

fit <- survfit(Surv(OS_m,Death) ~ group, data = cli)
pdf('/public/workspace/liuqzh/gastric_cancer/survival/GSE26253_State.pdf',height=5,width=5)
ggsurvplot(
	fit,
	data = cli,
	title = "Kaplan-Meier Survival Curve",
	pval = TRUE,  
	pval.method = TRUE,  
	pval.size = 6,  
	pval.coord = c(60, 0.1),  
	legend.title = "Group",
	xlab = "Time (Months)",
	ylab = "Survival Probability",
	font.main = 14,  
	font.x = 12,  
	font.y = 12,  
	font.tickslab = 10,  
	palette = c("#377EB8","#E41A1C"),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('GSE26253')
dev.off()
}

###########################
{
bytlib load gcc
bytlib load R-4.0.2
R

library(GSVA)
library(survival)
library(survminer)
library(dplyr)
 
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

Data=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv',header=T,sep=',')
Data_MP1<-list(Data$MP_1)
Data_MP2<-list(Data$MP_2)
Data_MP3<-list(Data$MP_3)

ssGSEA_Score_MP1=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP1, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP2=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP2, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP3=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP3, method='ssgsea', kcdf='Poisson')))#ssGSEA计算

colnames(ssGSEA_Score_MP1)='MP1'
colnames(ssGSEA_Score_MP2)='MP2'
colnames(ssGSEA_Score_MP3)='MP3'

ssGSEA_Score=cbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$MP1_Sub='Other'
ssGSEA_Score$MP1_Sub[which(ssGSEA_Score$MP1 > quantile(ssGSEA_Score$MP1,0.5))]='MP1'
ssGSEA_Score$MP2_Sub='Other'
ssGSEA_Score$MP2_Sub[which(ssGSEA_Score$MP2 > quantile(ssGSEA_Score$MP2,0.5))]='MP2'
ssGSEA_Score$MP3_Sub='Other'
ssGSEA_Score$MP3_Sub[which(ssGSEA_Score$MP3 > quantile(ssGSEA_Score$MP3,0.5))]='MP3'

ssGSEA_Score$Subtype_State='Variable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='MP1' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='MP2' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='MP3')]='Stable State'

ssGSEA_Score$sample=rownames(ssGSEA_Score)
cli$sample=cli$Sample

cli=merge(cli,ssGSEA_Score,by='sample')
rownames(cli)=cli$sample

cli$group<-factor(cli$Subtype_State,levels=c('Variable State','Stable State'))

fit <- survfit(Surv(OS,event) ~ group, data = cli)
pdf('/public/workspace/liuqzh/gastric_cancer/survival/GSE84426_State.pdf',height=5,width=5)
ggsurvplot(
	fit,
	data = cli,
	title = "Kaplan-Meier Survival Curve",
	pval = TRUE,  
	pval.method = TRUE,  
	pval.size = 6,  
	pval.coord = c(60, 0.1),  
	legend.title = "Group",
	xlab = "Time (Months)",
	ylab = "Survival Probability",
	font.main = 14,  
	font.x = 12,  
	font.y = 12,  
	font.tickslab = 10,  
	palette = c("#377EB8","#E41A1C"),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('GSE84426')
dev.off()
}

#############
bytlib load R-4.0.2
bytlib load gcc
R

library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)

GSE15459_data=readRDS('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE15459/GSE15459_data.rds')
GSE15459_info=read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE15459/GSE15459_data_clinical_information.csv',row.names=1)
GSE15459_info$sample=rownames(GSE15459_info)
cli=GSE15459_info

set.seed(12345)
dat_exp=GSE15459_data
dat_exp<-dat_exp[rowSums(dat_exp)>0,]

Data=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/MP_subtype_III.csv',header=T,sep=',')
Data_MP1<-list(Data$MP_1)
Data_MP2<-list(Data$MP_2)
Data_MP3<-list(Data$MP_3)

ssGSEA_Score_MP1=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP1, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP2=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP2, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score_MP3=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP3, method='ssgsea', kcdf='Poisson')))#ssGSEA计算

colnames(ssGSEA_Score_MP1)='MP1'
colnames(ssGSEA_Score_MP2)='MP2'
colnames(ssGSEA_Score_MP3)='MP3'

ssGSEA_Score=cbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2,ssGSEA_Score_MP3)
ssGSEA_Score=as.data.frame(ssGSEA_Score)

ssGSEA_Score$MP1_Sub='Other'
ssGSEA_Score$MP1_Sub[which(ssGSEA_Score$MP1 > quantile(ssGSEA_Score$MP1,0.5))]='MP1'
ssGSEA_Score$MP2_Sub='Other'
ssGSEA_Score$MP2_Sub[which(ssGSEA_Score$MP2 > quantile(ssGSEA_Score$MP2,0.5))]='MP2'
ssGSEA_Score$MP3_Sub='Other'
ssGSEA_Score$MP3_Sub[which(ssGSEA_Score$MP3 > quantile(ssGSEA_Score$MP3,0.5))]='MP3'

ssGSEA_Score$Subtype_State='Variable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='MP1' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='MP2' & ssGSEA_Score$MP3_Sub=='Other')]='Stable State'
ssGSEA_Score$Subtype_State[which(ssGSEA_Score$MP1_Sub=='Other' & ssGSEA_Score$MP2_Sub=='Other' & ssGSEA_Score$MP3_Sub=='MP3')]='Stable State'

ssGSEA_Score$sample=rownames(ssGSEA_Score)

cli=merge(cli,ssGSEA_Score,by='sample')
rownames(cli)=cli$sample

cli$group<-factor(cli$Subtype_State,levels=c('Variable State','Stable State'))

fit <- survfit(Surv(OS,Status) ~ group, data = cli)
pdf('/public/workspace/liuqzh/gastric_cancer/survival/GSE15459_State.pdf',height=5,width=5)
ggsurvplot(
	fit,
	data = cli,
	title = "Kaplan-Meier Survival Curve",
	pval = TRUE,  
	pval.method = TRUE,  
	pval.size = 6,  
	pval.coord = c(60, 0.1),  
	legend.title = "Group",
	xlab = "Time (Months)",
	ylab = "Survival Probability",
	font.main = 14,  
	font.x = 12,  
	font.y = 12,  
	font.tickslab = 10,  
	palette = c("#377EB8","#E41A1C"),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('GSE15459')
dev.off()



