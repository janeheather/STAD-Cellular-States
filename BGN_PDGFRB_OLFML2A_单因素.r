####OS####
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
Data_MP1<-list(c('BGN','PDGFRB','OLFML2A'))

dat_exp_3_marker=as.data.frame(t(dat_exp[which(rownames(dat_exp) %in% c('BGN','PDGFRB','OLFML2A')),]))
sur$sample=rownames(sur)
dat_exp_3_marker$sample=rownames(dat_exp_3_marker)

COX_data=merge(sur,dat_exp_3_marker,'sample')
rownames(COX_data)=COX_data$sample
COX_data=COX_data[,-c(1,3)]

df=COX_data
df$'OS.time'=as.numeric(df$'OS.time')
#设置p值的阈值
pfilter <- 0.05   
#新建空白数据框
uniresult <- data.frame()  
#使用for循环对输入数据中的100个基因依次进行单因素COX分析
#单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
for(i in colnames(df[,3:ncol(df)])){
	unicox <- coxph(Surv('OS.time','OS') ~ df[,i], data = df)
	unisum<- summary(unicox)   
	pvalue <- round(unisum$coefficients[,5],3) 
	if(pvalue<pfilter){ 
    uniresult <- rbind(uniresult,
                       cbind(gene=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
    ))
  }
}   


ssGSEA_Score=as.data.frame(t(gsva(as.matrix(dat_exp),Data_MP1, method='ssgsea', kcdf='Poisson')))#ssGSEA计算
ssGSEA_Score$sample=rownames(ssGSEA_Score)
colnames(ssGSEA_Score)=c('gene','sample')

sur$sample=rownames(sur)

cli=merge(sur,ssGSEA_Score,by='sample')
rownames(cli)=cli$sample
cli$group<-ifelse(cli$gene > median(cli$gene),"High","Low")
cli$group<-factor(cli$group,levels=c('Low','High'))

fit <- survfit(Surv(OS.time,OS) ~ group, data = cli)

pdf('/public/workspace/liuqzh/gastric_cancer/survival/TCGA_STAD_STr_score.pdf',height=5,width=5)
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


###############
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
Data=read.csv('/public/workspace/liuqzh/gastric_cancer/Endo_score/All_MP_subtype_Result.csv',header=T,sep=',')
genelist=c(Data$Endo_MP_5)
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
	
	HR=as.numeric(broom::tidy(model)[2])
	L95CI=as.numeric(broom::tidy(model)[3])
	H95CI=as.numeric(broom::tidy(model)[4])
	p_value=as.numeric(broom::tidy(model)[5])
	
	data_coxph=cbind(coeff,HR,L95CI,H95CI,p_value)
	
	result=rbind(result,data_coxph)
}
rownames(result)=colnames(sur_data[,-c(ncol(sur_data),c(ncol(sur_data)-1))])
result=as.data.frame(result)


