######
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)

########模块细胞定义#########

Cytokines_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Cytokines.txt',header=T,sep='\t')
Cytokines=Cytokines_list$Symbol

Antigen_Processing_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Antigen_Processing_and_Presentation.txt',header=T,sep='\t')
Antigen_Processing=Antigen_Processing_list$Symbol

Antimicrobials_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Antimicrobials.txt',header=T,sep='\t')
Antimicrobials=Antimicrobials_list$Symbol

Chemokine_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Chemokine_Receptors.txt',header=T,sep='\t')
Chemokine_Receptors=Chemokine_Receptors_list$Symbol

Cytokine_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Cytokine_Receptors.txt',header=T,sep='\t')
Cytokine_Receptors=Cytokine_Receptors_list$Symbol

Interferons_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Interferons.txt',header=T,sep='\t')
Interferons=Interferons_list$Symbol

Interferons_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Interferons_Receptors.txt',header=T,sep='\t')
Interferons_Receptors=Interferons_Receptors_list$Symbol

Natural_Killer_Cell_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Natural_Killer_Cell.txt',header=T,sep='\t')
Natural_Killer_Cell=Natural_Killer_Cell_list$Symbol

TGF_b_Family_Members_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/TGF_b_Family_Members.txt',header=T,sep='\t')
TGF_b_Family_Members=TGF_b_Family_Members_list$Symbol

TGF_b_Family_Members_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/TGF_b_Family_Members_Receptors.txt',header=T,sep='\t')
TGF_b_Family_Members_Receptors=TGF_b_Family_Members_Receptors_list$Symbol

TNF_Family_Members_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/TNF_Family_Members.txt',header=T,sep='\t')
TNF_Family_Members=TNF_Family_Members_list$Symbol

TNF_Family_Members_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/TNF_Family_Members_Receptors.txt',header=T,sep='\t')
TNF_Family_Members_Receptors=TNF_Family_Members_Receptors_list$Symbol

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

genelist=list(c('BGN','PDGFRB','OLFML2A'))
####
TPM<-TPM[rowSums(TPM)>0,]
GENE_list=list(TNF_Family_Members_Receptors)

ssGSEA_Score_MP1=t(gsva(as.matrix(TPM),genelist, method='ssgsea'))
ssGSEA_Score_MP2=t(gsva(as.matrix(TPM),GENE_list, method='ssgsea'))
ssGSEA_Score=cbind(ssGSEA_Score_MP1,ssGSEA_Score_MP2)
colnames(ssGSEA_Score)=c('Str','MP2')

ssGSEA_Score=as.data.frame(ssGSEA_Score)

data_f=ssGSEA_Score
p1=ggplot(data=data_f, aes(x=Str, y=MP2))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x=Str, y=MP2)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 1, label.y = 2.5) +  # 添加皮尔逊相关系数
	theme_minimal() +  # 使用简洁主题
	labs(
		title = "BGN and PDGFRB correlation",
		x = "Str",
		y = "MP2"
	) +  # 添加图形标题和轴标题
	theme(
		plot.title = element_text(hjust = 0.5),  # 标题居中
		axis.text = element_text(color = "gray20"),  # 轴文字颜色
		axis.title = element_text(color = "gray20"),  # 轴标题颜色
		panel.grid.major = element_blank(),  # 移除主要网格线
		panel.grid.minor = element_blank(),  # 移除次要网格线
		panel.background = element_rect(fill = "white", color = "gray50")  # 背景颜色
	)
p1
pdf('/public/workspace/liuqzh/gastric_cancer/figure3/BGN_and_PDGFRB_correlation.pdf',height=7,width=7)
p1
dev.off()






