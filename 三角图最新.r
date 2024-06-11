
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(dplyr)
library(ggplot2)
library(ggtern)

set.seed(520)
#读取scRNA
dat<-readRDS('/public/workspace/liuqzh/gastric/STAD_Epithelium_Integrate_celltype.rds')
DefaultAssay(dat) <- "RNA"
TPM=as.data.frame(dat[["RNA"]]@counts) %>% apply(2,function(x){x/sum(x) * 10000})
TPM<-TPM[rowSums(TPM)>0,]
TPM<-as.matrix(TPM)

genelist=read.csv("/public/workspace/liuqzh/gastric_cancer/scalop-master/subtype_III.csv",sep = ",",row.names = 1,header = F)
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

saveRDS(scdata,"/public/workspace/liuqzh/gastric_cancer/Module_Sample_Subtype_III.rds")
Module_Sample_Subtype_III<-readRDS("/public/workspace/liuqzh/gastric_cancer/Module_Sample_Subtype_III.rds")
data_meta=Module_Sample_Subtype_III

Subtype=NULL
for(i in 1:nrow(data_meta)){
	ch_max=pmax(data_meta[i,1],data_meta[i,2],data_meta[i,3])
	if(ch_max == data_meta[i,1]){
		Subtype=c(Subtype,'MP_1')
	}	
	if(ch_max == data_meta[i,2]){
		Subtype=c(Subtype,'MP_2')
	}	
	if(ch_max == data_meta[i,3]){
		Subtype=c(Subtype,'MP_3')
	}	
}
Subtype=as.matrix(Subtype)
rownames(Subtype)=rownames(Module_Sample_Subtype_III)
colnames(Subtype)='Subtype'

result=NULL
for(i in 1:nrow(data_meta)){
	EMT=c(0,data_meta[i,1])
	NLT=c(-data_meta[i,2]*cos(pi/6),-data_meta[i,2]*sin(pi/6))
	PLF=c(data_meta[i,3]*cos(pi/6),-data_meta[i,3]*sin(pi/6))
	P_N_E=EMT+NLT+PLF
	result=rbind(result,P_N_E)
}
result=as.matrix(result)

result[,1]=as.numeric(result[,1])
result[,2]=as.numeric(result[,2])

rownames(result)=rownames(data_meta)
colnames(result)=c("x_name","y_name")
result=as.data.frame(result)

library(ggplot2)
library(ggsci)

pdf('/public/workspace/liuqzh/gastric_cancer/Tumor/scrabble_score_Subtype_III_plot.pdf',width=9,height=8,useDingbats=F)
aix=result
aix=cbind(result,Subtype)

A <- c(-15, -10)
B <- c(15, -10)

AB_length <- sqrt((B[1] - A[1])^2 + (B[2] - A[2])^2)
height <- sqrt(AB_length^2 - (AB_length/2)^2)
C <- c((A[1] + B[1])/2, A[2] + height)

colors=c("MP_1"="#D1832C","MP_2"="#5B93D4","MP_3"="#A97598")
p1<-ggplot(aix, aes(x=x_name, y=y_name)) + 
	geom_point(aes(color = Subtype), size = 0.8) +
	stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", color = NA, alpha = 0.1) +# 添加密度趋势
	scale_fill_gradient(low = "#E5E1DC", high = "#BD3934") +
	scale_color_manual(values = colors)+
    theme(axis.ticks = element_blank(),
          panel.background = element_blank()) +
    ylab(NULL)+xlab(NULL)+ 
	geom_polygon(data = data.frame(x = c(A[1], B[1], C[1]), 
								   y = c(A[2], B[2], C[2]), 
                               group = c(1, 1, 1)),
               aes(x = x, y = y, group = group), 
               fill = NA, color = "black", linetype = "solid") 	+ 
	theme_minimal() +
	theme(axis.ticks = element_blank(), panel.grid = element_blank(), 
		axis.text = element_blank(), axis.title = element_blank()) +
    ylim(-18,18)+
    xlim(-18,18)	
	p1
dev.off()

