
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

dat=dat[,1:45]

pvalue<-apply(dat,1,function(x) wilcox.test(x[grep("^R",colnames(dat))],x[grep("^N",colnames(dat))],paired = F)$p.value)
log2FC<-apply(dat,1,function(x){
  log2( (mean(x[grep("^R",colnames(dat))])+1) / (mean(x[grep("^N",colnames(dat))])+1) ) 
})

dat$p<-pvalue
dat$log2FC<-log2FC

dat1<-filter(dat,pvalue<0.05,abs(log2FC)>1)
rownames(dat1[which(dat1$log2FC < -1),])

###
dat=dat[,-c(ncol(dat),(ncol(dat)-1))]
Data_MP1<-list(c('BGN','PDGFRB','OLFML2A'))
LY6E_neu_score=list(c('IFIT1','MX1','HERC5','IFI6','ISG15','IFIT3','RSAD2','GBP1','IFIT2','XAF1','PARP9','UBE2L6','IRF7','PARP14','APOL6'))
CAFs_markers=list(c('BGN','CHPF','IGFBP7','NTF3','HHIPL1','ST6GALNAC5'))

ssGSEA_Score=gsva(as.matrix(dat),Data_MP1, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)#ssGSEA计算
ssGSEA_Score=as.data.frame(ssGSEA_Score)
Result=as.data.frame(t(ssGSEA_Score))

library(ggplot2)
library(ggpubr)
library(reshape2)

Result$Group='Other'
Result$Group[grep('^N',rownames(Result))]='N'
Result$Group[grep('^R',rownames(Result))]='R'
colnames(Result)[1]='Score'

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
# 将数据从宽格式转换为长格式
df=Result
df$Group <- factor(df$Group, levels = c("N", "R"))  # 设置分组的顺序

df=as.data.frame(df)

p <- ggplot(df, aes(x = Group, y = Score, color = Group)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("red" = "red", "lightblue" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "PRJEB25780", x = "Group", y = "Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("N","R")), method='wilcox.test',
                            label = "p.signif")

pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/PRJEB25780/PRJEB25780_data_N_R.pdf',height=6,width=3)
print(p)
dev.off()


#ssGSEA_Score=ssGSEA_Score[,1:45,drop=F]
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#ssGSEA_Score<-normalize(ssGSEA_Score)
####得分差异
pvalue<-apply(ssGSEA_Score,1,function(x) wilcox.test(x[grep("^R",colnames(dat))],x[grep("^N",colnames(dat))],paired = F)$p.value)
pvalue
####
ssGSEA_Score=as.data.frame(t(ssGSEA_Score))
colnames(ssGSEA_Score)='Score'

ssGSEA_Score$group=''
ssGSEA_Score$group[grep('^N_',rownames(ssGSEA_Score))]='N'
ssGSEA_Score$group[grep('^R_',rownames(ssGSEA_Score))]='R'

library(pROC)
pdf('/public/workspace/liuqzh/gastric_cancer/bulk_data/PRJEB25780/PRJEB25780_data_N_R_ROC.pdf',height=4,width=4)
roc(ssGSEA_Score$group, ssGSEA_Score$Score, levels=c("R","N"), direction="<", plot=T, print.auc=TRUE, print.thres=TRUE)
dev.off()

##############


