bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)

Maligant_MP_1<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_1.txt',header=T,row.names=1)
Maligant_MP_2<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_2.txt',header=T,row.names=1)
Maligant_MP_3<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_3.txt',header=T,row.names=1)

normalize <- function(x) {
return((x - min(x)) / (max(x) - min(x)))
}
Maligant_MP_1=arrange(Maligant_MP_1,desc(Simplicity))
Maligant_MP_1$Simplicity_score=normalize(Maligant_MP_1$Simplicity)
Maligant_MP_2=arrange(Maligant_MP_2,desc(Simplicity))
Maligant_MP_2$Simplicity_score=normalize(Maligant_MP_2$Simplicity)
Maligant_MP_3=arrange(Maligant_MP_3,desc(Simplicity))
Maligant_MP_3$Simplicity_score=normalize(Maligant_MP_3$Simplicity)

colnames(Maligant_MP_1)[15]='ADDS'
colnames(Maligant_MP_2)[15]='ADDS'
colnames(Maligant_MP_3)[15]='ADDS'

DATA_MP=rbind(Maligant_MP_1,Maligant_MP_2,Maligant_MP_3)
DATA_MP$MP1 = -log10(DATA_MP$Maligant_MP_1_pval)
DATA_MP$MP2 = -log10(DATA_MP$Maligant_MP_2_pval)
DATA_MP$MP3 = -log10(DATA_MP$Maligant_MP_3_pval)

DATA_MP$simplicity='Low'
DATA_MP$simplicity[which(DATA_MP$Simplicity_score >= 0.8)]='High'

tcga_est<-read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/TCGA/tcga_est.csv',header=T,row.names=1)
tcga_est=as.data.frame(t(tcga_est))
tcga_est=tcga_est[-1,]
rownames(tcga_est)=gsub('\\.','-',rownames(tcga_est))

DATA_MP$sample=rownames(DATA_MP)
tcga_est$sample=rownames(tcga_est)

DATA_MP_EST=merge(DATA_MP,tcga_est,'sample')
rownames(DATA_MP_EST)=DATA_MP_EST$sample
DATA_MP_EST=as.data.frame(DATA_MP_EST)
DATA_MP_EST=DATA_MP_EST[,-1]
DATA_MP_EST$Subtype=factor(DATA_MP_EST$Subtype,levels=c('MP_1','MP_2','MP_3'))
DATA_MP_EST$Stromal_Score=as.numeric(DATA_MP_EST$StromalScore)
DATA_MP_EST$Immune_Score=as.numeric(DATA_MP_EST$ImmuneScore)
DATA_MP_EST$ESTIMATE_Score=as.numeric(DATA_MP_EST$ESTIMATEScore)
DATA_MP_EST$Tumor_Purity_Score=as.numeric(DATA_MP_EST$TumorPurity)
#####
p <- ggplot(DATA_MP_EST, aes(x = Subtype, y = Stromal_Score, color = Subtype)) +
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
	labs(title = "TCGA-STAD", x = "Group", y = "Stromal_Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c('MP_1','MP_2'),c('MP_2','MP_3'),c('MP_1','MP_3')), method='wilcox.test',
                            label = "p.signif")
pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_Stromal_Score.pdf",width=4,height=6)
print(p)
dev.off()
#####
###
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)

Maligant_MP_1<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_1.txt',header=T,row.names=1)
Maligant_MP_2<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_2.txt',header=T,row.names=1)
Maligant_MP_3<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_3.txt',header=T,row.names=1)

normalize <- function(x) {
return((x - min(x)) / (max(x) - min(x)))
}
Maligant_MP_1=arrange(Maligant_MP_1,desc(Simplicity))
Maligant_MP_1$Simplicity_score=normalize(Maligant_MP_1$Simplicity)
Maligant_MP_2=arrange(Maligant_MP_2,desc(Simplicity))
Maligant_MP_2$Simplicity_score=normalize(Maligant_MP_2$Simplicity)
Maligant_MP_3=arrange(Maligant_MP_3,desc(Simplicity))
Maligant_MP_3$Simplicity_score=normalize(Maligant_MP_3$Simplicity)

colnames(Maligant_MP_1)[15]='ADDS'
colnames(Maligant_MP_2)[15]='ADDS'
colnames(Maligant_MP_3)[15]='ADDS'

DATA_MP=rbind(Maligant_MP_1,Maligant_MP_2,Maligant_MP_3)
DATA_MP$MP1 = -log10(DATA_MP$Maligant_MP_1_pval)
DATA_MP$MP2 = -log10(DATA_MP$Maligant_MP_2_pval)
DATA_MP$MP3 = -log10(DATA_MP$Maligant_MP_3_pval)

DATA_MP$simplicity='Low'
DATA_MP$simplicity[which(DATA_MP$Simplicity_score >= 0.8)]='High'

tcga_est<-read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/TCGA/tcga_est.csv',header=T,row.names=1)
tcga_est=as.data.frame(t(tcga_est))
tcga_est=tcga_est[-1,]
rownames(tcga_est)=gsub('\\.','-',rownames(tcga_est))

DATA_MP$sample=rownames(DATA_MP)
tcga_est$sample=rownames(tcga_est)

DATA_MP_EST=merge(DATA_MP,tcga_est,'sample')
rownames(DATA_MP_EST)=DATA_MP_EST$sample
DATA_MP_EST=as.data.frame(DATA_MP_EST)
DATA_MP_EST=DATA_MP_EST[,-1]
DATA_MP_EST$Subtype=factor(DATA_MP_EST$Subtype,levels=c('MP_1','MP_2','MP_3'))
DATA_MP_EST$Stromal_Score=as.numeric(DATA_MP_EST$StromalScore)
DATA_MP_EST$Immune_Score=as.numeric(DATA_MP_EST$ImmuneScore)
DATA_MP_EST$ESTIMATE_Score=as.numeric(DATA_MP_EST$ESTIMATEScore)
DATA_MP_EST$Tumor_Purity_Score=as.numeric(DATA_MP_EST$TumorPurity)

DATA_MP_EST$Tumor_Purity='High'
DATA_MP_EST$Tumor_Purity[which(DATA_MP_EST$Tumor_Purity_Score <= 0.7)]='Low'
DATA_MP_EST$PATIENT=rownames(DATA_MP_EST)

sur<-read.table('/public/workspace/yumiao/STAD/tcga/TCGA-STAD.survival.tsv',header=T,sep='\t',row.names=1)
sur$OS.time=sur$OS.time/30
sur$OS.time[which(sur$OS.time >= 60)]=60
sur$OS[which(sur$OS.time >= 50)]=0
sur$PATIENT=paste(sur$PATIENT,'-01A',sep='')

DATA_MP_sur=merge(sur,DATA_MP_EST,'PATIENT')
fit <- survfit(Surv(OS.time,OS) ~ Tumor_Purity, data = DATA_MP_sur)

pdf('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_sur.pdf',height=5,width=5)
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
dev.off()


###
#####
p <- ggplot(DATA_MP_EST, aes(x = Subtype, y = Immune_Score, color = Subtype)) +
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
	labs(title = "TCGA-STAD", x = "Group", y = "Immune_Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c('MP_1','MP_2'),c('MP_2','MP_3'),c('MP_1','MP_3')), method='wilcox.test',
                            label = "p.signif")
pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_Immune_Score.pdf",width=4,height=6)
print(p)
dev.off()
#####

#####
p <- ggplot(DATA_MP_EST, aes(x = Subtype, y = ESTIMATE_Score, color = Subtype)) +
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
	labs(title = "TCGA-STAD", x = "Group", y = "ESTIMATE_Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c('MP_1','MP_2'),c('MP_2','MP_3'),c('MP_1','MP_3')), method='wilcox.test',
                            label = "p.signif")
pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_ESTIMATE_Score.pdf",width=4,height=6)
print(p)
dev.off()
#####

#####
p <- ggplot(DATA_MP_EST, aes(x = Subtype, y = Tumor_Purity_Score, color = Subtype)) +
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
	labs(title = "TCGA-STAD", x = "Group", y = "Tumor_Purity_Score")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c('MP_1','MP_2'),c('MP_2','MP_3'),c('MP_1','MP_3')), method='wilcox.test',
                            label = "p.signif")
pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_Tumor_Purity_Score.pdf",width=4,height=6)
print(p)
dev.off()
#####

bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)

Maligant_MP_1<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_1.txt',header=T,row.names=1)
Maligant_MP_2<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_2.txt',header=T,row.names=1)
Maligant_MP_3<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_3.txt',header=T,row.names=1)

normalize <- function(x) {
return((x - min(x)) / (max(x) - min(x)))
}
Maligant_MP_1=arrange(Maligant_MP_1,desc(Simplicity))
Maligant_MP_1$Simplicity_score=normalize(Maligant_MP_1$Simplicity)
Maligant_MP_2=arrange(Maligant_MP_2,desc(Simplicity))
Maligant_MP_2$Simplicity_score=normalize(Maligant_MP_2$Simplicity)
Maligant_MP_3=arrange(Maligant_MP_3,desc(Simplicity))
Maligant_MP_3$Simplicity_score=normalize(Maligant_MP_3$Simplicity)

colnames(Maligant_MP_1)[15]='ADDS'
colnames(Maligant_MP_2)[15]='ADDS'
colnames(Maligant_MP_3)[15]='ADDS'

DATA_MP=rbind(Maligant_MP_1,Maligant_MP_2,Maligant_MP_3)
DATA_MP$MP1 = -log10(DATA_MP$Maligant_MP_1_pval)
DATA_MP$MP2 = -log10(DATA_MP$Maligant_MP_2_pval)
DATA_MP$MP3 = -log10(DATA_MP$Maligant_MP_3_pval)

DATA_MP$simplicity='Low'
DATA_MP$simplicity[which(DATA_MP$Simplicity_score >= 0.8)]='High'
DATA_MP$PATIENT=rownames(DATA_MP)

sur<-read.table('/public/workspace/yumiao/STAD/tcga/TCGA-STAD.survival.tsv',header=T,sep='\t',row.names=1)
sur$OS.time=sur$OS.time/30
sur$OS.time[which(sur$OS.time >= 60)]=60
sur$OS[which(sur$OS.time >= 50)]=0
sur$PATIENT=paste(sur$PATIENT,'-01A',sep='')

DATA_MP_sur=merge(sur,DATA_MP,'PATIENT')
#DATA_MP_sur=DATA_MP_sur[which(DATA_MP_sur$simplicity == 'Low'),]
###

fit <- survfit(Surv(OS.time,OS) ~ Subtype, data = DATA_MP_sur)

pdf('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_sur.pdf',height=5,width=5)
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
dev.off()

#######
#DATA_MP_sur=DATA_MP_sur[which(DATA_MP_sur$simplicity == 'High'),]
DATA_MP_sur$non_MP3='MP3'
DATA_MP_sur$non_MP3[which(DATA_MP_sur$Subtype %in% c('MP_1','MP_2'))]='non_MP3'
fit <- survfit(Surv(OS.time,OS) ~ non_MP3, data = DATA_MP_sur)
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

write.table(DATA_MP,"/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt",row.names=T,col.names=T,quote=F,sep="\t")






