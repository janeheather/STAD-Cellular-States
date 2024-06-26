
setwd('/public/workspace/liuqzh/gastric_cancer/figure3/')
library(CellChat)
library(Seurat)
library(future)
cellchat<-readRDS('/public/workspace/yumiao/STAD/NMF_class_cellchat/cellchat.rds')
pdf('inter2.pdf',height=10)
p = netVisual_bubble(cellchat,remove.isolate = T,sources.use = c('Maligant_MP_2'),targets.use =c('Endo_MP_5'),signaling = c('PDGF','EGF','TGFb','GALECTIN'))
print(p)
dev.off()