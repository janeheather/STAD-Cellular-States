########三角图
###t2_5931
bytlib load R-4.0.2
bytlib load gcc
R

library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)
library(Seurat)
library(dplyr)
library(ggtern)
Module_Sample_Subtype_III<-readRDS("/public/workspace/liuqzh/gastric_cancer/Module_Sample_Subtype_III.rds")
data_meta=Module_Sample_Subtype_III
data_meta<-data_meta[grep('^t2_5931',rownames(data_meta)),]

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
rownames(Subtype)=rownames(data_meta)
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

pdf('/public/workspace/liuqzh/gastric_cancer/Tumor/scrabble_score_Subtype_III_plot_t2_5931.pdf',width=9,height=8,useDingbats=F)
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


###P17_EGC
bytlib load R-4.0.2
bytlib load gcc
R

library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)
library(Seurat)
library(dplyr)
library(ggtern)
Module_Sample_Subtype_III<-readRDS("/public/workspace/liuqzh/gastric_cancer/Module_Sample_Subtype_III.rds")
data_meta=Module_Sample_Subtype_III
data_meta<-data_meta[grep('^P17_EGC',rownames(data_meta)),]

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
rownames(Subtype)=rownames(data_meta)
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

pdf('/public/workspace/liuqzh/gastric_cancer/Tumor/scrabble_score_Subtype_III_plot_P17_EGC.pdf',width=9,height=8,useDingbats=F)
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

###Tu13
bytlib load R-4.0.2
bytlib load gcc
R

library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)
library(Seurat)
library(dplyr)
library(ggtern)
Module_Sample_Subtype_III<-readRDS("/public/workspace/liuqzh/gastric_cancer/Module_Sample_Subtype_III.rds")
data_meta=Module_Sample_Subtype_III
data_meta<-data_meta[grep('^Tu13',rownames(data_meta)),]

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
rownames(Subtype)=rownames(data_meta)
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

pdf('/public/workspace/liuqzh/gastric_cancer/Tumor/scrabble_score_Subtype_III_plot_Tu13.pdf',width=9,height=8,useDingbats=F)
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

###Tu15
bytlib load R-4.0.2
bytlib load gcc
R

library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggpubr)
library(rlang)
library(Seurat)
library(dplyr)
library(ggtern)
Module_Sample_Subtype_III<-readRDS("/public/workspace/liuqzh/gastric_cancer/Module_Sample_Subtype_III.rds")
data_meta=Module_Sample_Subtype_III
data_meta<-data_meta[grep('^Tu15',rownames(data_meta)),]

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
rownames(Subtype)=rownames(data_meta)
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

pdf('/public/workspace/liuqzh/gastric_cancer/Tumor/scrabble_score_Subtype_III_plot_Tu15.pdf',width=9,height=8,useDingbats=F)
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




