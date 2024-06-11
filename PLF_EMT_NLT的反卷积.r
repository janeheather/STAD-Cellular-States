#########PLF_EMT_NLT的反卷积###########
bytlib load R-4.0.2
bytlib load gcc
R

library(e1071)
library(preprocessCore)
library(parallel)
library(CIBERSORT)
library(Seurat)

STAD_FPKM=readRDS("/public/workspace/liuqzh/gastric/STAD_FPKM_data_exp.rds")
genelist<-read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/subtype_III.csv",sep = ",",row.names = 1,header = F)
rownames(genelist)=c('PLF','EMT','NLT')

EMT <-as.character(genelist[2,])
EMT<-EMT[EMT != ""]
NLT=as.character(genelist[3,])
NLT<-NLT[NLT != ""]
PLF <-as.character(genelist[1,])
PLF<-PLF[PLF != ""]

dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Module_CNA_Sample_step_III_Subtype_III.rds')
Module_Sample_Subtype_III<-readRDS("/public/workspace/liuqzh/gastric_cancer/Module_Sample_Subtype_III.rds")

Module_Sample_Subtype_III=Module_Sample_Subtype_III[intersect(colnames(dat),rownames(Module_Sample_Subtype_III)),]
dat$PLF=Module_Sample_Subtype_III$MP_1
dat$EMT=Module_Sample_Subtype_III$MP_2
dat$NLT=Module_Sample_Subtype_III$MP_3

Idents(dat)<-dat$Subtype
averger_exp <- AverageExpression(dat)[[1]]
averger_exp=averger_exp[which(as.data.frame(apply(averger_exp,1,mean))!=0),]

all_genes=c(EMT,PLF,NLT)
averger=averger_exp[which(rownames(averger_exp) %in% all_genes),]

sig_matrix=averger

test_expr=as.matrix(STAD_FPKM)
#write.table(test_expr, "/public/workspace/liuqzh/gastric/test_expr_results.txt", row.names = T, sep = "\t")

results <- cibersort(sig_matrix = sig_matrix, mixture_file = test_expr)
write.table(results, "/public/workspace/liuqzh/gastric/EMT_PLF_NLT_results_all.txt", row.names = T, sep = "\t")


rowscale <- results[,1:ncol(sig_matrix)]#只是相当于备份了一下results
rowscale <- rowscale[,apply(rowscale, 2, function(x){sum(x)>0})]#删除全是0的列
pheatmap(rowscale,
         scale = 'row',#按行标准化，不标准化就会按绝对值显示，很诡异
         cluster_col=T,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
         cluster_row=T,#是否对行聚类
         angle_col = "315")#调整X轴坐标的倾斜角度