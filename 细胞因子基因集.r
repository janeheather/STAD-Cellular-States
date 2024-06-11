######
bytlib load gcc
bytlib load R-4.0.2
R

library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)

########模块细胞定义#########

Cytokines_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Cytokines.txt',header=T,sep='\t')
Cytokines=Cytokines_list$Symbol

Antigen_Processing_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Antigen_Processing_and_Presentation.txt',header=T,sep='\t')
Antigen_Processing=Antigen_Processing_list$Symbol

Antimicrobials_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Antimicrobials.txt',header=T,sep='\t')
Antimicrobials=Antimicrobials_list$Symbol

Chemokine_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Chemokine_Receptors.txt',header=T,sep='\t')
Chemokine_Receptors=Chemokine_Receptors_list$Symbol

Cytokine_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Cytokine_Receptors.txt',header=T,sep='\t')
Cytokine_Receptors=Cytokine_Receptors_list$Symbol

Interferons_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Interferons.txt',header=T,sep='\t')
Interferons=Interferons_list$Symbol

Interferons_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Interferons_Receptors.txt',header=T,sep='\t')
Interferons_Receptors=Interferons_Receptors_list$Symbol

Natural_Killer_Cell_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/Natural_Killer_Cell.txt',header=T,sep='\t')
Natural_Killer_Cell=Natural_Killer_Cell_list$Symbol

TGF_b_Family_Members_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/TGF_b_Family_Members.txt',header=T,sep='\t')
TGF_b_Family_Members=TGF_b_Family_Members_list$Symbol

TGF_b_Family_Members_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/TGF_b_Family_Members_Receptors.txt',header=T,sep='\t')
TGF_b_Family_Members_Receptors=TGF_b_Family_Members_Receptors_list$Symbol

TNF_Family_Members_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/TNF_Family_Members.txt',header=T,sep='\t')
TNF_Family_Members=TNF_Family_Members_list$Symbol

TNF_Family_Members_Receptors_list=read.table('/public/workspace/liuqzh/gastric_cancer/immport_genelist/TNF_Family_Members_Receptors.txt',header=T,sep='\t')
TNF_Family_Members_Receptors=TNF_Family_Members_Receptors_list$Symbol









