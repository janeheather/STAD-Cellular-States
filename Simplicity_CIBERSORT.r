bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)

CIBERSORT_Results<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/CIBERSORT_Results_LM22.csv',header=T,sep=',',row.names=1)
colnames(CIBERSORT_Results)=gsub('\\.','_',colnames(CIBERSORT_Results))
rownames(CIBERSORT_Results)=paste(rownames(CIBERSORT_Results),'-01A',sep='')
CIBERSORT_Results$PATIENT=rownames(CIBERSORT_Results)

DATA_MP<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t',row.names=1)

CIBERSORT_Results_DATA_MP=merge(CIBERSORT_Results,DATA_MP,'PATIENT')
rownames(CIBERSORT_Results_DATA_MP)=CIBERSORT_Results_DATA_MP$PATIENT
CIBERSORT_Results_DATA_MP=CIBERSORT_Results_DATA_MP[,-1]

CIBERSORT=CIBERSORT_Results_DATA_MP[,c(1:22,38)]
CIBERSORT_plot<-melt(CIBERSORT)

pdf('/public/workspace/liuqzh/gastric_cancer/Simplicity/CIBERSORT_plot.pdf',width=12,height=3,useDingbats=F)
ggplot(CIBERSORT_plot, aes(x=variable,y=value,colour=Subtype)) +
    stat_boxplot(geom="errorbar")+
    geom_boxplot(outlier.shape = NA, notch=FALSE)+
#   geom_point(shape=21, position = position_jitter(0.25), size=1.5, aes(fill=Subtype), colour="black")+
	scale_colour_manual(values=c('MP_1'='#D1832C','MP_2'='#5B92D3','MP_3'='#A97597'))+
    stat_compare_means(aes(group = Subtype), method = "t.test",label = "p.signif",hide.ns = TRUE)+
	theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1), panel.border=element_rect(color='black', fill=NA), panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"), legend.position="right")+
    labs(x="")
dev.off()
