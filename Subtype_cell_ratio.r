bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
#读取scRNA
dat<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Epithelium_Subtype_Score.rds')
DefaultAssay(dat) <- "RNA"

data_ratio=as.data.frame(table(dat$seurat_clusters,dat$Subtype))
colnames(data_ratio)[1:2]=c('Clusters','Subtype')
ratio=data_ratio %>% group_by(Clusters) %>% mutate(percent = Freq/sum(Freq)) %>% as.data.frame()
ratio$Clusters<-paste0('C',ratio$Clusters)
ratio<-ratio[,c(1,2,4)]
c<-spread(ratio,key=Subtype,value=percent)
c$a<-gsub('C','',c$Clusters) %>% as.numeric
c<-arrange(c,a)
rownames(c)<-c[,1]
c<-c[,-c(1,5)]


cell_marker<-FindAllMarkers(dat,logfc.threshold=0.5,only.pos = TRUE, min.pct = 0.25)
top10 <- cell_marker%>% group_by(cluster) %>% top_n(10, avg_log2FC)
a<-NULL
for(i in 0:28){
	tmp_gene<-list(filter(top10,cluster==i)$gene)
	tmp_dat<-subset(dat,subset=seurat_clusters==i)
	tmp_score<-AddModuleScore(tmp_dat,features = tmp_gene)
	b<-c(cor(tmp_score$MP_1,tmp_score$Cluster1),cor.test(tmp_score$MP_1,tmp_score$Cluster1)$p.value,cor(tmp_score$MP_2,tmp_score$Cluster1),cor.test(tmp_score$MP_2,tmp_score$Cluster1)$p.value,cor(tmp_score$MP_3,tmp_score$Cluster1),cor.test(tmp_score$MP_3,tmp_score$Cluster1)$p.value)

	a<-rbind(a,b)
}
a<-as.data.frame(a)
rownames(a)<-paste0('C',0:28)
colnames(a)<-c('MP_1_COR','MP_1_P','MP_2_COR','MP_2_P','MP_3_COR','MP_3_P')

aaa<-a[,c(1,3,5)]
aaa$id<-rownames(aaa)
aaa$id<-factor(aaa$id,levels=paste0("C",0:28))
aaa<-melt(aaa)
c$id<-rownames(c)
d<-melt(c)
d$id<-factor(d$id,levels=paste0("C",0:28))

aaa$size<-d$value
p = ggplot(aaa, aes(id,variable)) +
	geom_tile(aes(fill = value),colour = "white") + 
	scale_fill_gradient2(mid='white',low = "darkblue",high = "red")+
	geom_point(aes(size=size),color="black",shape=21,alpha=1)

pdf('/public/workspace/liuqzh/gastric_cancer/ste_score/subtype_cell_ratio.pdf',width=12,height=3,useDingbats=F)
print(p)
dev.off()









