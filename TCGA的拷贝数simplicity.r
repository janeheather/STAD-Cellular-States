##########################################################################################
#########PLF_EMT_NLT的反卷积###########
bytlib load R-4.0.2
bytlib load gcc
R

library(e1071)
library(preprocessCore)
library(parallel)
library(CIBERSORT)
library(Seurat)
library(TCGAbiolinks)
library(dplyr)

DATA_MP=read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/DATA_MP.txt',header=T,sep='\t')

DATA_MP$simplicity='Low'
DATA_MP$simplicity[which(DATA_MP$Simplicity_score >= 0.8)]='High'
DATA_MP=DATA_MP[which(DATA_MP$simplicity == 'High'),]

setwd('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/')
query <- GDCquery(project = "TCGA-STAD", 
                  data.category = "Copy Number Variation", 
                  data.type = "Masked Copy Number Segment")

GDCdownload(query, method = "api", files.per.chunk = 100)
segment_dat <- GDCprepare(query = query)
segment_dat$Sample <- substring(segment_dat$Sample,1,16)
segment_dat <- grep("01A$",segment_dat$Sample) %>% segment_dat[.,]
segment_dat[,1] <- segment_dat$Sample
segment_dat <- segment_dat[,-7]
segment_dat<-as.data.frame(segment_dat)

MP_1<-segment_dat[which(segment_dat[,1] %in% rownames(DATA_MP[which(DATA_MP$Subtype == 'MP_1'),])),]
MP_2<-segment_dat[which(segment_dat[,1] %in% rownames(DATA_MP[which(DATA_MP$Subtype == 'MP_2'),])),]
MP_3<-segment_dat[which(segment_dat[,1] %in% rownames(DATA_MP[which(DATA_MP$Subtype == 'MP_3'),])),]

write.table(MP_1,"/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_1/MaskedCopyNumberSegment_MP_1.txt",sep="\t",
            quote = F,col.names = F,row.names = F)
write.table(MP_2,"/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_2/MaskedCopyNumberSegment_MP_2.txt",sep="\t",
            quote = F,col.names = F,row.names = F)		
write.table(MP_3,"/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_3/MaskedCopyNumberSegment_MP_3.txt",sep="\t",
            quote = F,col.names = F,row.names = F)

snp<-read.table('/public/workspace/yumiao/STAD/snp6.na35.remap.hg38.subset.txt',header=T)
snp<-snp[snp$freqcnv=='FALSE',]
snp<-snp[,1:3]
colnames(snp)<-c("Marker name","chromosome","Marker position")
write.table(snp,'/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/Marker.txt',sep="\t",
            quote = F,col.names = F,row.names = F)

#linux

basedir=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_1/
segfile=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_1/MaskedCopyNumberSegment_MP_1.txt
markersfile=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/Marker.txt
refgenefile=/public/workspace/yumiao/STAD/gistic2_result/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

/public/workspace/yumiao/STAD/gistic2_result/gistic2 -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.95 -armpeel 1 -savegene 1 -gcm extreme
##########

basedir=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_2/
segfile=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_2/MaskedCopyNumberSegment_MP_2.txt
markersfile=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/Marker.txt
refgenefile=/public/workspace/yumiao/STAD/gistic2_result/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

/public/workspace/yumiao/STAD/gistic2_result/gistic2 -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.95 -armpeel 1 -savegene 1 -gcm extreme
##########

basedir=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_3/
segfile=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_3/MaskedCopyNumberSegment_MP_3.txt
markersfile=/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/Marker.txt
refgenefile=/public/workspace/yumiao/STAD/gistic2_result/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

/public/workspace/yumiao/STAD/gistic2_result/gistic2 -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.95 -armpeel 1 -savegene 1 -gcm extreme
##########



R
setwd('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_1/')
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)#BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
                 )
df$chromNum <- 1:length(df$chromName) # 染色体名字纯数字版，为了和scores.gistic里面的对应

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息

df

df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df
scores <- read.table("scores.gistic",sep="\t",header=T,stringsAsFactors = F)

head(scores)
chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
range(scores$G.score)
## [1] 0.000000 1.038373

scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1

range(scores$G.score)
## [1] -0.564793  0.280607
library(ggplot2)
library(ggsci)

## 使用maftools分析
library(maftools)

scores$'G.score'=scores$'X.log10.q.value.'*scores$'G.score'

pdf("/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_1/MP_1_cnv.pdf",width=12,height=4)
df$ypos <- rep(c(0.2,0.25),11)
ggplot(scores, aes(StartPos, G.score))+
  geom_line(aes(group=Type, fill=factor(Type,levels = c("Del","Amp")),colour = Type))+
  scale_color_manual(values = c('#BA3630','#5B92D3'))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-5,5)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
        )
dev.off()

########
########

R
setwd('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_2/')
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)#BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
                 )
df$chromNum <- 1:length(df$chromName) # 染色体名字纯数字版，为了和scores.gistic里面的对应

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息

df

df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df
scores <- read.table("scores.gistic",sep="\t",header=T,stringsAsFactors = F)

head(scores)
chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
range(scores$G.score)
## [1] 0.000000 1.038373

scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1

range(scores$G.score)
## [1] -0.564793  0.280607
library(ggplot2)
library(ggsci)

## 使用maftools分析
library(maftools)

scores$'G.score'=scores$'X.log10.q.value.'*scores$'G.score'

pdf("/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_2/MP_2_cnv.pdf",width=12,height=4)
df$ypos <- rep(c(0.2,0.25),11)
ggplot(scores, aes(StartPos, G.score))+
  geom_line(aes(group=Type, fill=factor(Type,levels = c("Del","Amp")),colour = Type))+
  scale_color_manual(values = c('#BA3630','#5B92D3'))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-5,5)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
        )
dev.off()

##########
##########

R
setwd('/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_3/')
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)#BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
                 )
df$chromNum <- 1:length(df$chromName) # 染色体名字纯数字版，为了和scores.gistic里面的对应

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息

df

df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df
scores <- read.table("scores.gistic",sep="\t",header=T,stringsAsFactors = F)

head(scores)
chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
range(scores$G.score)
## [1] 0.000000 1.038373

scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1

range(scores$G.score)
## [1] -0.564793  0.280607
library(ggplot2)
library(ggsci)

## 使用maftools分析
library(maftools)

scores$'G.score'=scores$'X.log10.q.value.'*scores$'G.score'

pdf("/public/workspace/liuqzh/gastric_cancer/CNV/simplicity/MP_3/MP_3_cnv.pdf",width=12,height=4)
df$ypos <- rep(c(0.2,0.25),11)
ggplot(scores, aes(StartPos, G.score))+
  geom_line(aes(group=Type, fill=factor(Type,levels = c("Del","Amp")),colour = Type))+
  scale_color_manual(values = c('#BA3630','#5B92D3'))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-5,5)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
        )
dev.off()


