bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(reshape2)
library(dplyr)
library(pheatmap)
library(philentropy)
library(tibble)
library(tidyr)
library(patchwork)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(stringr)
library("Rtsne")
library(RColorBrewer)
library(scales)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ggplot2)
library(ggpubr)
library(rlang)
source('/public/workspace/liuqzh/ssgseaMOD.r')

setwd('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP1')
####
sur<-read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_YUHS_sur.csv',header=T,row.names=2)
gse<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_matrix.txt',header=T,sep='\t',row.names=1)
ref<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GPL6884.txt',header=T,sep='\t',row.names=1)

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

marker<-read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv',header=T,row.names=1)
All_Gene=c(marker$Maligant_MP_1,marker$Maligant_MP_2,marker$Maligant_MP_3)

###
Gene=marker$Maligant_MP_1

exp=dat_exp[which(rownames(dat_exp) %in% All_Gene),]

mod.generate(Gene,out=paste0('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP1/','Maligant_MP_1','.mod'))
mod <- mod.analyze2(exp,'Maligant_MP_1','/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP1/',permN=100000)
rownames(mod)=gsub("\\.","-",rownames(mod))

write.table(mod,file="/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP1/Maligant_MP_1_GSE62254.txt",row.names=T,col.names=T,quote=F,sep="\t")
##

########
########
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(reshape2)
library(dplyr)
library(pheatmap)
library(philentropy)
library(tibble)
library(tidyr)
library(patchwork)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(stringr)
library("Rtsne")
library(RColorBrewer)
library(scales)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ggplot2)
library(ggpubr)
library(rlang)
source('/public/workspace/liuqzh/ssgseaMOD.r')

setwd('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP2')
####
sur<-read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_YUHS_sur.csv',header=T,row.names=2)
gse<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_matrix.txt',header=T,sep='\t',row.names=1)
ref<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GPL6884.txt',header=T,sep='\t',row.names=1)

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

marker<-read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv',header=T,row.names=1)
All_Gene=c(marker$Maligant_MP_1,marker$Maligant_MP_2,marker$Maligant_MP_3)

###
Gene=marker$Maligant_MP_2

exp=dat_exp[which(rownames(dat_exp) %in% All_Gene),]

mod.generate(Gene,out=paste0('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP2/','Maligant_MP_2','.mod'))
mod <- mod.analyze2(exp,'Maligant_MP_2','/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP2/',permN=100000)
rownames(mod)=gsub("\\.","-",rownames(mod))

write.table(mod,file="/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP2/Maligant_MP_2_GSE62254.txt",row.names=T,col.names=T,quote=F,sep="\t")
##

########
########
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(reshape2)
library(dplyr)
library(pheatmap)
library(philentropy)
library(tibble)
library(tidyr)
library(patchwork)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(stringr)
library("Rtsne")
library(RColorBrewer)
library(scales)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)
library(ggplot2)
library(ggpubr)
library(rlang)
source('/public/workspace/liuqzh/ssgseaMOD.r')

setwd('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP3')
####
sur<-read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_YUHS_sur.csv',header=T,row.names=2)
gse<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_matrix.txt',header=T,sep='\t',row.names=1)
ref<-read.table('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GPL6884.txt',header=T,sep='\t',row.names=1)

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

marker<-read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv',header=T,row.names=1)
All_Gene=c(marker$Maligant_MP_1,marker$Maligant_MP_2,marker$Maligant_MP_3)

###
Gene=marker$Maligant_MP_3

exp=dat_exp[which(rownames(dat_exp) %in% All_Gene),]

mod.generate(Gene,out=paste0('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP3/','Maligant_MP_3','.mod'))
mod <- mod.analyze2(exp,'Maligant_MP_3','/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP3/',permN=100000)
rownames(mod)=gsub("\\.","-",rownames(mod))

write.table(mod,file="/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP3/Maligant_MP_3_GSE62254.txt",row.names=T,col.names=T,quote=F,sep="\t")
##

bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(reshape2)
library(dplyr)
library(pheatmap)
library(philentropy)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(stringr)
library("Rtsne")
library(RColorBrewer)
library(scales)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)

Maligant_MP_1<-read.table('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP1/Maligant_MP_1_GSE13861.txt',header=T,row.names=1)
Maligant_MP_2<-read.table('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP2/Maligant_MP_2_GSE13861.txt',header=T,row.names=1)
Maligant_MP_3<-read.table('/public/workspace/liuqzh/ssgsea_tmp/GSE13861/MP3/Maligant_MP_3_GSE13861.txt',header=T,row.names=1)
Maligant_MP_1$Rank_MP1=rank(Maligant_MP_1$Maligant_MP_1_norm)
Maligant_MP_2$Rank_MP2=rank(Maligant_MP_2$Maligant_MP_2_norm)
Maligant_MP_3$Rank_MP3=rank(Maligant_MP_3$Maligant_MP_3_norm)
TCGA_subtype=cbind(Maligant_MP_1,Maligant_MP_2,Maligant_MP_3)

data_meta=TCGA_subtype
# 定义函数用于找到最小的pval和相应的排名
find_min_pval <- function(Maligant_MP_1_pval, Maligant_MP_2_pval, Maligant_MP_3_pval, Rank_MP1, Rank_MP2, Rank_MP3) {
  pvals <- c(Maligant_MP_1_pval, Maligant_MP_2_pval, Maligant_MP_3_pval)
  ranks <- c(Rank_MP1, Rank_MP2, Rank_MP3)
  
  min_pval_indices <- which(pvals == min(pvals))
  
  if (length(min_pval_indices) == 1) {
    return(min_pval_indices)
  } else {
    return(min_pval_indices[which.max(ranks[min_pval_indices])])
  }
}

# 使用pmap函数找到每行的MP类型
data_meta$Subtype <- pmap_chr(data_meta, function(Maligant_MP_1_pval, Maligant_MP_2_pval, Maligant_MP_3_pval, Rank_MP1, Rank_MP2, Rank_MP3, ...) {
  min_index <- find_min_pval(Maligant_MP_1_pval, Maligant_MP_2_pval, Maligant_MP_3_pval, Rank_MP1, Rank_MP2, Rank_MP3)
  return(c("MP_1", "MP_2", "MP_3")[min_index])
})

data_meta=as.data.frame(data_meta)	

Maligant_MP_1=data_meta[which(data_meta$Subtype == 'MP_1'),]
Maligant_MP_2=data_meta[which(data_meta$Subtype == 'MP_2'),]
Maligant_MP_3=data_meta[which(data_meta$Subtype == 'MP_3'),]

#####MP1
ADNS=NULL
for(i in rownames(data_meta[which(data_meta$Subtype == 'MP_1'),])){
	Rank_MP2=data_meta[which(rownames(data_meta) == i),'Rank_MP2']
	Pval_MP2=data_meta[which(rownames(data_meta) == i),'Maligant_MP_2_pval']	
	MP2_tmp=data_meta[which(data_meta$Rank_MP2 < Rank_MP2),]
	ADNS_1=sum(MP2_tmp$Maligant_MP_2_pval)-(Rank_MP2-1)*Pval_MP2

	Rank_MP3=data_meta[which(rownames(data_meta) == i),'Rank_MP3']
	Pval_MP3=data_meta[which(rownames(data_meta) == i),'Maligant_MP_3_pval']	
	MP3_tmp=data_meta[which(data_meta$Rank_MP3 < Rank_MP3),]
	ADNS_2=sum(MP3_tmp$Maligant_MP_3_pval)-(Rank_MP3-1)*Pval_MP3
	
	ADNS_ALL=ADNS_1+ADNS_2
	ADNS=c(ADNS,ADNS_ALL)
}
Maligant_MP_1$ADNS=ADNS
#######
R0_MP_1=min(data_meta$Maligant_MP_1_pval)
Rn1_MP_1=max(data_meta$Maligant_MP_1_pval)
MP_1_ADDS=sum(data_meta$Maligant_MP_1_pval)+nrow(data_meta)*R0_MP_1
Maligant_MP_1=cbind(Maligant_MP_1,MP_1_ADDS)
Maligant_MP_1$Simplicity=(Maligant_MP_1$MP_1_ADDS-Maligant_MP_1$ADNS)*(Rn1_MP_1-R0_MP_1)/nrow(Maligant_MP_1)
write.table(Maligant_MP_1, "/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_1_GSE13861.txt", row.names = T, sep = "\t")


#####MP2
ADNS=NULL
for(i in rownames(data_meta[which(data_meta$Subtype == 'MP_2'),])){
	Rank_MP1=data_meta[which(rownames(data_meta) == i),'Rank_MP1']
	Pval_MP1=data_meta[which(rownames(data_meta) == i),'Maligant_MP_1_pval']	
	MP1_tmp=data_meta[which(data_meta$Rank_MP1 < Rank_MP1),]
	ADNS_1=sum(MP1_tmp$Maligant_MP_1_pval)-(Rank_MP1-1)*Pval_MP1

	Rank_MP3=data_meta[which(rownames(data_meta) == i),'Rank_MP3']
	Pval_MP3=data_meta[which(rownames(data_meta) == i),'Maligant_MP_3_pval']	
	MP3_tmp=data_meta[which(data_meta$Rank_MP3 < Rank_MP3),]
	ADNS_2=sum(MP3_tmp$Maligant_MP_3_pval)-(Rank_MP3-1)*Pval_MP3
	
	ADNS_ALL=ADNS_1+ADNS_2
	ADNS=c(ADNS,ADNS_ALL)
}
Maligant_MP_2$ADNS=ADNS
#######
R0_MP_2=min(data_meta$Maligant_MP_2_pval)
Rn1_MP_2=max(data_meta$Maligant_MP_2_pval)
MP_2_ADDS=sum(data_meta$Maligant_MP_2_pval)+nrow(data_meta)*R0_MP_2
Maligant_MP_2=cbind(Maligant_MP_2,MP_2_ADDS)
Maligant_MP_2$Simplicity=(Maligant_MP_2$MP_2_ADDS-Maligant_MP_2$ADNS)*(Rn1_MP_2-R0_MP_2)/nrow(Maligant_MP_2)
write.table(Maligant_MP_2, "/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_2_GSE13861.txt", row.names = T, sep = "\t")


#####MP3
ADNS=NULL
for(i in rownames(data_meta[which(data_meta$Subtype == 'MP_3'),])){
	Rank_MP1=data_meta[which(rownames(data_meta) == i),'Rank_MP1']
	Pval_MP1=data_meta[which(rownames(data_meta) == i),'Maligant_MP_1_pval']	
	MP1_tmp=data_meta[which(data_meta$Rank_MP1 < Rank_MP1),]
	ADNS_1=sum(MP1_tmp$Maligant_MP_1_pval)-(Rank_MP1-1)*Pval_MP1

	Rank_MP2=data_meta[which(rownames(data_meta) == i),'Rank_MP2']
	Pval_MP2=data_meta[which(rownames(data_meta) == i),'Maligant_MP_2_pval']	
	MP2_tmp=data_meta[which(data_meta$Rank_MP2 < Rank_MP2),]
	ADNS_2=sum(MP2_tmp$Maligant_MP_2_pval)-(Rank_MP2-1)*Pval_MP2
	
	ADNS_ALL=ADNS_1+ADNS_2
	ADNS=c(ADNS,ADNS_ALL)
}
Maligant_MP_3$ADNS=ADNS
#######
R0_MP_3=min(data_meta$Maligant_MP_3_pval)
Rn1_MP_3=max(data_meta$Maligant_MP_3_pval)
MP_3_ADDS=sum(data_meta$Maligant_MP_3_pval)+nrow(data_meta)*R0_MP_3
Maligant_MP_3=cbind(Maligant_MP_3,MP_3_ADDS)
Maligant_MP_3$Simplicity=(Maligant_MP_3$MP_3_ADDS-Maligant_MP_3$ADNS)*(Rn1_MP_3-R0_MP_3)/nrow(Maligant_MP_3)
write.table(Maligant_MP_3, "/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_3_GSE13861.txt", row.names = T, sep = "\t")

#################
############
module load gcc
module load R-4.0.2
R


library(Seurat)
library(scater)
library(reshape2)
library(dplyr)
library(pheatmap)
library(philentropy)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(survival)
library(survminer)
library(GSVA)
library(estimate)
library(corrplot)
library(ComplexHeatmap)
library(stringr)
library("Rtsne")
library(RColorBrewer)
library(scales)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)

normalize <- function(x) {
return((x - min(x)) / (max(x) - min(x)))
}

Maligant_MP_1<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_1_GSE13861.txt',header=T,row.names=1)
Maligant_MP_2<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_2_GSE13861.txt',header=T,row.names=1)
Maligant_MP_3<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_3_GSE13861.txt',header=T,row.names=1)
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

DATA_MP$GEO_ID=rownames(DATA_MP)
sur<-read.csv('/public/workspace/liuqzh/gastric_cancer/bulk_data/GSE13861/GSE13861_YUHS_sur.csv',header=T,row.names=2)
sur$GEO_ID=rownames(sur)

sur=sur[intersect(rownames(DATA_MP),rownames(sur)),]
DATA_MP=cbind(DATA_MP,sur)

#write.table(DATA_MP,file='/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_GSE13861.txt',row.names=T,col.names=T,quote=F,sep="\t")

pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_GSE13861.pdf",width=12,height=3)
library(circlize)
DATA_MP$simplicity='Low'
DATA_MP$simplicity[which(DATA_MP$Simplicity_score >= 0.8)]='High'

ha = HeatmapAnnotation(Subtype = as.character(DATA_MP$'Subtype'),
					   simplicity = as.character(DATA_MP$'simplicity'),
					   Subgroup = as.character(DATA_MP$'Subgroup'),
					   Location = as.character(DATA_MP$'Location'),
					   Lauren = as.character(DATA_MP$'Lauren'),
					   AJCC6 = as.character(DATA_MP$'AJCC6'),
					   Adjuvant_chemo = as.character(DATA_MP$'Adjuvant_chemo'),
					   col = list(simplicity = c("Low" = "#427D39",'High'="#BA3630"),
								  Subtype= c("MP_1"='#D1832C',"MP_2"='#5B92D3',"MP_3"='#A97597'),
								  Subgroup = c('EP' = '#CF812A', 'MP' = '#5A92D2'),
								  Location = c('antrum' = '#BA3630','body'='#CF812A','cardia'='#5A92D2','entire'='#D1832C','fundus'='#A97597'),
								  Lauren = c('diffuse' = '#9F70A6','intestinal' = '#427D39','mixed'='#282177'),
								  AJCC6 = c("I" = "#BA3630", "II" = "#9F70A6",'III'="#427D39",'IV'='#70174F'),
								  Adjuvant_chemo = c('Yes' = '#BA3630','No' = '#427D39')))

mycols <- colorRamp2(breaks = c(0,1.5,5),
                    colors = c('#5A8FCA','#F2F2F0','#E31A1C'))
Heatmap(t(DATA_MP[,c(18:20)]),cluster_columns = F,cluster_rows = F,show_column_names = F,name = "heat",top_annotation = ha,
		col = mycols)
dev.off()

cli=DATA_MP
fit <- survfit(Surv(OS_m,Death) ~ Subtype, data = cli)
summary(fit)$table

#cli$non_MP3='MP3'
#cli$non_MP3[which(cli$Subtype %in% c('MP_1','MP_2'))]='non_MP3'
#fit <- survfit(Surv(OS_m,Death) ~ non_MP3, data = cli)

#cli$non_MP3='MP1'
#cli$non_MP3[which(cli$Subtype %in% c('MP_2','MP_3'))]='non_MP1'
#fit <- survfit(Surv(OS_m,Death) ~ non_MP3, data = cli)

#cli$non_MP3='MP2'
#cli$non_MP3[which(cli$Subtype %in% c('MP_1','MP_3'))]='non_MP2'
#fit <- survfit(Surv(OS_m,Death) ~ non_MP3, data = cli)

pdf('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_GSE13861_sur.pdf',height=5,width=5)
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
	palette = c('#D1832C','#5B92D3','#A97597'),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('GSE13861')
dev.off()

#####
fit <- survfit(Surv(RFS_m,Recurrence) ~ Subtype, data = cli)
summary(fit)$table

#cli$non_MP3='MP3'
#cli$non_MP3[which(cli$Subtype %in% c('MP_1','MP_2'))]='non_MP3'
#fit <- survfit(Surv(RFS_m,Recurrence) ~ non_MP3, data = cli)

#cli$non_MP3='MP1'
#cli$non_MP3[which(cli$Subtype %in% c('MP_2','MP_3'))]='non_MP1'
#fit <- survfit(Surv(RFS_m,Recurrence) ~ non_MP3, data = cli)

#cli$non_MP3='MP2'
#cli$non_MP3[which(cli$Subtype %in% c('MP_1','MP_3'))]='non_MP2'
#fit <- survfit(Surv(RFS_m,Recurrence) ~ non_MP3, data = cli)

pdf('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_GSE13861_RFS.pdf',height=5,width=5)
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
	palette = c('#D1832C','#5B92D3','#A97597'),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('GSE13861')
dev.off()

pdf("/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP_Simplicity_GSE13861.pdf",width=12,height=3)
Simplicity_colors <- colorRampPalette(c("#ABA19A", "#E3E0DD"))(100)
Heatmap(t(DATA_MP[,17]),col = Simplicity_colors,cluster_columns = F,cluster_rows = F,show_column_names = F)
dev.off()


cli$non_MP3='MP3'
cli$non_MP3[which(cli$Subtype %in% c('MP_1','MP_2'))]='non_MP3'
fit <- survfit(Surv(OS_m,Death) ~ non_MP3, data = cli)

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
	palette = c('#D1832C','#5B92D3','#A97597'),  
	censor.shape = 124,  
	censor.size = 3  
)+ggtitle('GSE13861')


