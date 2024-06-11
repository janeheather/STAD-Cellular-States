bytlib load R-4.0.2
bytlib load gcc
R

library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

mycds<-readRDS('/public/workspace/liuqzh/gastric_cancer/E_A/simple_p/mycds_STAD_P_E_T_N_cell_SeuratFeatures.rds')
monocle_meta_data=mycds@phenoData@data
Cho_Epi=read.csv("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/Cho_Epithelium.auc_mtx.csv",sep=",",header=T,row.names=1)
colnames(Cho_Epi)=gsub("\\...","",colnames(Cho_Epi))

Cho_Epi=Cho_Epi[which(rownames(Cho_Epi) %in% rownames(monocle_meta_data)),]
monocle_meta_data=monocle_meta_data[which(rownames(Cho_Epi) %in% rownames(monocle_meta_data)),]

STAD_P_RNA_Cho_Epi=cbind(Cho_Epi,monocle_meta_data$Pseudotime)
colnames(STAD_P_RNA_Cho_Epi)[309]='Pseudotime'

###PLF
STAD_P_RNA_Cho_Epi=STAD_P_RNA_Cho_Epi[which(rownames(STAD_P_RNA_Cho_Epi) %in% rownames(monocle_meta_data[which(monocle_meta_data$subtype %in% c('NLT','TSL','pre-PLF','PLF')),])),]

M_Cho_Epi=cor(STAD_P_RNA_Cho_Epi,method='pearson')
Pseudotime_Cho_Epi=M_Cho_Epi[which(colnames(M_Cho_Epi) == 'Pseudotime'),-309]
Pseudotime_Cho_Epi_list=names(which(Pseudotime_Cho_Epi >= 0.05))

[1] "CEBPG" "FOXF1" "NR1I2" "NR1I3"

###EMT
STAD_P_RNA_Cho_Epi=STAD_P_RNA_Cho_Epi[which(rownames(STAD_P_RNA_Cho_Epi) %in% rownames(monocle_meta_data[which(monocle_meta_data$subtype %in% c('NLT','TSL','pre-EMT','EMT')),])),]

M_Cho_Epi=cor(STAD_P_RNA_Cho_Epi,method='pearson')
Pseudotime_Cho_Epi=M_Cho_Epi[which(colnames(M_Cho_Epi) == 'Pseudotime'),-309]
Pseudotime_Cho_Epi_list=names(which(Pseudotime_Cho_Epi >= 0.05))

 [1] "ARID3A"  "BACH1"   "CDX1"    "CDX2"    "CEBPG"   "CREB3L2" "CUX1"
 [8] "ELF4"    "ESRRA"   "ETS2"    "ETV3"    "FOSL1"   "FOSL2"   "FOXO3"
[15] "HNF4A"   "HNF4G"   "HOXB13"  "IRF7"    "JUND"    "KLF16"   "KLF3"
[22] "KLF5"    "KLF6"    "MAF"     "MAFG"    "MAFK"    "MLXIPL"  "MXD1"
[29] "NFE2L1"  "NR1H4"   "NR1I2"   "NR1I3"   "NR2F6"   "NR5A2"   "PPARA"
[36] "PPARD"   "PRDM1"   "RREB1"   "RXRA"    "SMAD3"   "SP1"     "SP2"
[43] "STAT2"   "TBX10"   "TBX3"    "TCF7L2"  "TFE3"    "THRB"    "VDR"
[50] "ZBTB7B"  "ZNF165"  "ZNF91"

###NLT
STAD_P_RNA_Cho_Epi=STAD_P_RNA_Cho_Epi[which(rownames(STAD_P_RNA_Cho_Epi) %in% rownames(monocle_meta_data[which(monocle_meta_data$subtype %in% c('NLT','TSL','pre-PLF','PLF','pre-EMT','EMT')),])),]

M_Cho_Epi=cor(STAD_P_RNA_Cho_Epi,method='pearson')
Pseudotime_Cho_Epi=M_Cho_Epi[which(colnames(M_Cho_Epi) == 'Pseudotime'),-309]
Pseudotime_Cho_Epi_list=names(which(Pseudotime_Cho_Epi <= -0.05))

[1] "BCL11A" "E2F8"   "GATA6"  "NFIC"   "NR2F2"  "SPIB"   "TEAD2"
