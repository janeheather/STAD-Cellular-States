###P13_EGC
{
bytlib load R-4.0.2
R

library(Seurat)
library(velocyto.R)
library(tidyverse)
library(SeuratWrappers)
##数据基础分析
##数据基础分析
# 读取loom文件

STAD_P <- readRDS("/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/STAD_P_E_T_N.rds")
DefaultAssay(STAD_P)<-'RNA'

STAD_P=subset(STAD_P,cells=which(STAD_P$subtype=='PLF'|STAD_P$subtype=='pre-PLF'|
								 STAD_P$subtype=='EMT'|STAD_P$subtype=='pre-EMT'|
								 STAD_P$subtype=='TSL'|STAD_P$subtype=='NLT'))

DAT=subset(STAD_P,cells=which(STAD_P$clas=='P01_EGC'))
velo <- readRDS("/public/workspace/liuqzh/gastric_cancer/E_A/Epithelium/velocyto/STAD_Epithelium_MS_loom.rds")

PNET=DAT
cells=intersect(colnames(velo),colnames(PNET))

velo=velo[,cells]
PNET=PNET[,cells]

Idents(velo)=PNET$subtype
velo$celltype=PNET$subtype
velo@graphs=PNET@graphs
velo@reductions=PNET@reductions

##给细胞分配颜色
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = velo)))
ident.colors <- c("#70655C","#377EB8","#984EA3") #####修改颜色
names(x = ident.colors) <- levels(x = velo)
cell.colors <- ident.colors[Idents(object = velo)]
names(x = cell.colors) <- colnames(x = velo)

##速率分析
velo <- RunVelocity(velo, deltaT = 1, kCells = 25, fit.quantile = 0.02, 
        spliced.average = 0.2, unspliced.average = 0.05, ncores = 18)
#kCells：用于斜率平滑度计算最近邻细胞数量，越大越平滑，越小越能反映局部特征
#fit.quantile：0.02代表对基因表达量最高2%与最低2%的值执行gamma拟合
#spliced.average：过滤低表达丰度基因的标准，计算的是基因在cluster内的平均counts值
#unspliced.average：同上

##全局速率可视化
emb = Embeddings(PNET, reduction = "umap")
vel = Tool(velo, slot = "RunVelocity")

setwd("/public/workspace/liuqzh/gastric_cancer/E_A/Epithelium/")
pdf('STAD_P13_EGC.pdf')
show.velocity.on.embedding.cor(emb = emb, vel = vel, n = 200, scale = "sqrt", 
       cell.colors = ac(cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, 
       show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
       arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
dev.off()
}