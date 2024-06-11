bytlib load R-4.0.2
bytlib load gcc
R

library(parallel)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rlang)

DATA_MP=read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t')
DATA_MP$subgroup=paste(DATA_MP$simplicity,DATA_MP$Subtype,sep='_')

dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

TCGA_exp=dat_exp[,intersect(colnames(dat_exp),rownames(DATA_MP))]
TCGA_exp=as.matrix(TCGA_exp)
TCGA_exp_fibmarker=as.data.frame(t(TCGA_exp[which(rownames(TCGA_exp) %in% c('COL1A1','COL1A2','COL3A1','COL5A2','COL6A1','COL6A2','SPARC')),]))

TCGA_exp_fibmarker$PATIENT=rownames(TCGA_exp_fibmarker)
DATA_MP_plot=merge(DATA_MP,TCGA_exp_fibmarker,'PATIENT')

rownames(DATA_MP_plot)=DATA_MP_plot$PATIENT
DATA_MP_plot$subgroup=factor(DATA_MP_plot$subgroup,levels=c('High_MP_1','Low_MP_1','High_MP_2','Low_MP_2','High_MP_3','Low_MP_3'))

df=DATA_MP_plot

for(i in c('COL1A1','COL1A2','COL3A1','COL5A2','COL6A1','COL6A2','SPARC')){
y_names=as.name(i)
p <-ggplot(df, aes(x = subgroup, y = {{y_names}}, color = subgroup)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("High" = "red", "Low" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "TCGA-STAD", x = "Group", y = paste("Expression of",{{y_names}},sep=' '))
	p <- p + stat_compare_means(comparisons = list(c('High_MP_3','High_MP_1'),c('High_MP_3','Low_MP_1'),c('High_MP_3','High_MP_2'),c('High_MP_3','Low_MP_2'),c('High_MP_3','Low_MP_3')), method='wilcox.test',label = "p.signif")

pdf(paste('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/',{{y_names}},'_exp.pdf',sep=''),height=6,width=5)
print(p)
dev.off()
}

##########
bytlib load R-4.0.2
bytlib load gcc
R

library(parallel)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rlang)

DATA_MP=read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t')
DATA_MP$subgroup=paste(DATA_MP$simplicity,DATA_MP$Subtype,sep='_')

dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

TCGA_exp=dat_exp[,intersect(colnames(dat_exp),rownames(DATA_MP))]
TCGA_exp=as.matrix(TCGA_exp)
TCGA_exp_fibmarker=as.data.frame(t(TCGA_exp[which(rownames(TCGA_exp) %in% c('CTLA4','HAVCR2','LAG3','PDCD1','TIGIT')),]))

TCGA_exp_fibmarker$PATIENT=rownames(TCGA_exp_fibmarker)
DATA_MP_plot=merge(DATA_MP,TCGA_exp_fibmarker,'PATIENT')

rownames(DATA_MP_plot)=DATA_MP_plot$PATIENT
DATA_MP_plot$subgroup=factor(DATA_MP_plot$subgroup,levels=c('High_MP_1','Low_MP_1','High_MP_2','Low_MP_2','High_MP_3','Low_MP_3'))

df=DATA_MP_plot

for(i in c('CTLA4','HAVCR2','LAG3','PDCD1','TIGIT')){
y_names=as.name(i)
p <-ggplot(df, aes(x = subgroup, y = {{y_names}}, color = subgroup)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("High" = "red", "Low" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "TCGA-STAD", x = "Group", y = paste("Expression of",{{y_names}},sep=' '))
	p <- p + stat_compare_means(comparisons = list(c('High_MP_1','Low_MP_1'),c('High_MP_1','High_MP_2'),c('High_MP_1','Low_MP_2'),c('High_MP_1','High_MP_3'),c('High_MP_1','Low_MP_3')), method='wilcox.test',label = "p.signif")

pdf(paste('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/',{{y_names}},'_exp.pdf',sep=''),height=6,width=5)
print(p)
dev.off()
}

##########
bytlib load R-4.0.2
bytlib load gcc
R

library(parallel)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rlang)

DATA_MP=read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t')
DATA_MP$subgroup=paste(DATA_MP$simplicity,DATA_MP$Subtype,sep='_')

dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

TCGA_exp=dat_exp[,intersect(colnames(dat_exp),rownames(DATA_MP))]
TCGA_exp=as.matrix(TCGA_exp)
TCGA_exp_fibmarker=as.data.frame(t(TCGA_exp[which(rownames(TCGA_exp) %in% c('CTLA4','HAVCR2','LAG3','PDCD1','TIGIT')),]))

TCGA_exp_fibmarker$PATIENT=rownames(TCGA_exp_fibmarker)
DATA_MP_plot=merge(DATA_MP,TCGA_exp_fibmarker,'PATIENT')

rownames(DATA_MP_plot)=DATA_MP_plot$PATIENT
DATA_MP_plot$subgroup=factor(DATA_MP_plot$subgroup,levels=c('High_MP_1','Low_MP_1','High_MP_2','Low_MP_2','High_MP_3','Low_MP_3'))

df=DATA_MP_plot

for(i in c('CTLA4','HAVCR2','LAG3','PDCD1','TIGIT')){
y_names=as.name(i)
p <-ggplot(df, aes(x = subgroup, y = {{y_names}}, color = subgroup)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("High" = "red", "Low" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "TCGA-STAD", x = "Group", y = paste("Expression of",{{y_names}},sep=' '))
	p <- p + stat_compare_means(comparisons = list(c('High_MP_1','Low_MP_1'),c('High_MP_1','High_MP_2'),c('High_MP_1','Low_MP_2'),c('High_MP_1','High_MP_3'),c('High_MP_1','Low_MP_3')), method='wilcox.test',label = "p.signif")

pdf(paste('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/',{{y_names}},'_exp.pdf',sep=''),height=6,width=5)
print(p)
dev.off()
}

#############
#############
##########
bytlib load R-4.0.2
bytlib load gcc
R

library(parallel)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rlang)

DATA_MP=read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t')
DATA_MP$subgroup=paste(DATA_MP$simplicity,DATA_MP$Subtype,sep='_')

dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

TCGA_exp=dat_exp[,intersect(colnames(dat_exp),rownames(DATA_MP))]
TCGA_exp=as.matrix(TCGA_exp)
TCGA_exp_fibmarker=as.data.frame(t(TCGA_exp[which(rownames(TCGA_exp) %in% c('CXCL1','CXCL2','CXCL3','CXCL5','CCL20',
																			'SLC12A2','CCL3','CCL4','CXCL6','CXCL8',
																			'CXCL9','CXCL10','CXCL11')),]))

TCGA_exp_fibmarker$PATIENT=rownames(TCGA_exp_fibmarker)
DATA_MP_plot=merge(DATA_MP,TCGA_exp_fibmarker,'PATIENT')

rownames(DATA_MP_plot)=DATA_MP_plot$PATIENT
DATA_MP_plot$subgroup=factor(DATA_MP_plot$subgroup,levels=c('High_MP_1','Low_MP_1','High_MP_2','Low_MP_2','High_MP_3','Low_MP_3'))

df=DATA_MP_plot

for(i in c('CXCL1','CXCL2','CXCL3','CXCL5','CCL20',
		   'SLC12A2','CCL3','CCL4','CXCL6','CXCL8',
		   'CXCL9','CXCL10','CXCL11')){
y_names=as.name(i)
p <-ggplot(df, aes(x = subgroup, y = {{y_names}}, color = subgroup)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("High" = "red", "Low" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "TCGA-STAD", x = "Group", y = paste("Expression of",{{y_names}},sep=' '))
	p <- p + stat_compare_means(comparisons = list(c('High_MP_1','Low_MP_1'),c('High_MP_2','Low_MP_2'),c('High_MP_3','Low_MP_3')), method='wilcox.test',label = "p.signif")

pdf(paste('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/',{{y_names}},'_exp.pdf',sep=''),height=6,width=5)
print(p)
dev.off()
}

##########
bytlib load R-4.0.2
bytlib load gcc
R

library(parallel)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rlang)

DATA_MP=read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t')
DATA_MP$subgroup=paste(DATA_MP$simplicity,DATA_MP$Subtype,sep='_')

dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

TCGA_exp=dat_exp[,intersect(colnames(dat_exp),rownames(DATA_MP))]
TCGA_exp=as.matrix(TCGA_exp)
TCGA_exp_fibmarker=as.data.frame(t(TCGA_exp[which(rownames(TCGA_exp) %in% c('AQP5','LGR5')),]))

TCGA_exp_fibmarker$PATIENT=rownames(TCGA_exp_fibmarker)
DATA_MP_plot=merge(DATA_MP,TCGA_exp_fibmarker,'PATIENT')

rownames(DATA_MP_plot)=DATA_MP_plot$PATIENT
DATA_MP_plot$subgroup=factor(DATA_MP_plot$subgroup,levels=c('High_MP_1','Low_MP_1','High_MP_2','Low_MP_2','High_MP_3','Low_MP_3'))

df=DATA_MP_plot

for(i in c('AQP5','LGR5')){
y_names=as.name(i)
p <-ggplot(df, aes(x = subgroup, y = {{y_names}}, color = subgroup)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("High" = "red", "Low" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "TCGA-STAD", x = "Group", y = paste("Expression of",{{y_names}},sep=' '))
	p <- p + stat_compare_means(comparisons = list(c('High_MP_1','Low_MP_1'),c('High_MP_1','High_MP_2'),c('High_MP_1','Low_MP_2'),c('High_MP_1','High_MP_3'),c('High_MP_1','Low_MP_3')), method='wilcox.test',label = "p.signif")

pdf(paste('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/',{{y_names}},'_exp.pdf',sep=''),height=6,width=5)
print(p)
dev.off()
}

##########
bytlib load R-4.0.2
bytlib load gcc
R

library(parallel)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(rlang)

DATA_MP=read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t')
DATA_MP$subgroup=paste(DATA_MP$simplicity,DATA_MP$Subtype,sep='_')

dat_exp<-read.csv('/public/workspace/yumiao/STAD/tcga/stad_tcga_fpkm.csv',header=T,row.names=1)
colnames(dat_exp)=gsub('\\.','-',colnames(dat_exp))
dat_exp=dat_exp[,grep("-01A",colnames(dat_exp))]

TCGA_exp=dat_exp[,intersect(colnames(dat_exp),rownames(DATA_MP))]
TCGA_exp=as.matrix(TCGA_exp)
TCGA_exp_fibmarker=as.data.frame(t(TCGA_exp[which(rownames(TCGA_exp) %in% c('AQP5','CXCL1','CXCL3','IL1B','IL6')),]))

TCGA_exp_fibmarker$PATIENT=rownames(TCGA_exp_fibmarker)
DATA_MP_plot=merge(DATA_MP,TCGA_exp_fibmarker,'PATIENT')

rownames(DATA_MP_plot)=DATA_MP_plot$PATIENT
DATA_MP_plot$simplicity=factor(DATA_MP_plot$simplicity,levels=c('High','Low'))

df=DATA_MP_plot

for(i in c('AQP5','CXCL1','CXCL3','IL1B','IL6')){
y_names=as.name(i)
p <-ggplot(df, aes(x = simplicity, y = {{y_names}}, color = simplicity)) +
	geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
	geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
	scale_color_manual(values = c("High" = "red", "Low" = "lightblue")) +  # 设置颜色
	theme_minimal() +  # 使用简洁主题
	theme(  # 调整主题去掉不必要的元素
		panel.grid.major = element_blank(),  # 移除网格线
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),  # 移除背景
		axis.line = element_line(color = "black"),  # 坐标轴线
		legend.position = "none"  # 不显示图例
	) +
	labs(title = "TCGA-STAD", x = "Group", y = paste("Expression of",{{y_names}},sep=' '))
	p <- p + stat_compare_means(comparisons = list(c('High','Low')), method='wilcox.test',label = "p.signif")

pdf(paste('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/','High_Low_',{{y_names}},'_exp.pdf',sep=''),height=5,width=3)
print(p)
dev.off()
}


#################
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

DATA_MP=read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t')
DATA_MP$sur_subtype='NULL'
DATA_MP$sur_subtype[which(DATA_MP$Subtype %in% c('MP_1','MP_2'))]='non_MP_3'
DATA_MP$sur_subtype[which(DATA_MP$Subtype %in% c('MP_3'))]='MP_3'

sur<-read.table('/public/workspace/yumiao/STAD/tcga/TCGA-STAD.survival.tsv',header=T,sep='\t',row.names=1)
sur$OS.time=sur$OS.time/30
sur$OS.time[which(sur$OS.time >= 60)]=60
sur$OS[which(sur$OS.time >= 50)]=0
sur$PATIENT=paste(sur$PATIENT,'-01A',sep='')

DATA_MP_sur=merge(sur,DATA_MP,'PATIENT')
fit <- survfit(Surv(OS.time,OS) ~ sur_subtype, data = DATA_MP_sur)

ggsurvplot(
	fit,
	data = DATA_MP_sur,
	title = "Kaplan-Meier Survival Curve",
	pval = TRUE,  
	pval.method = TRUE,  
	pval.size = 6,  
	pval.coord = c(30, 0.1),  
	legend.title = "Group",
	xlab = "Time (Months)",
	ylab = "Survival Probability",
	font.main = 14,  
	font.x = 12,  
	font.y = 12,  
	font.tickslab = 10,  
	palette = c('#D1832C','#5B92D3','#A97597'),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('TCGA_STAD')



