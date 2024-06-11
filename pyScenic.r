module load hdf5-1.8.13
module load r/4.0.4
module load gcc
R
.libPaths("/public/workspace/liuqzh/R/library/4.0.4/")
##############################################
# run pyscenic
##############################################
options(stringsAsFactors = FALSE)
library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(stringr)
library(readr)
library(SeuratDisk)
library(dplyr)
library(tidyr)
seurat2loom <- function(inRDS, outloom){
 dat.loom <- SeuratDisk::as.loom(inRDS, filename = outloom, assay = "RNA", verbose = FALSE, overwrite = TRUE)
    dat.loom$close_all()
}

STAD_P<-readRDS('/public/workspace/liuqzh/gastric_cancer/Intergrate_data/STAD_Integrate_celltype.rds')
DefaultAssay(STAD_P)<-'RNA'
{
setwd("/public/workspace/liuqzh/gastric_cancer/pySCENIC/fib/")
#########################fibroblast
Cho_Epi=subset(STAD_P,cells=which(STAD_P$celltype=='Fibroblast'))
DefaultAssay(Cho_Epi) <- "RNA"

Cho_Epi <- DietSeurat(Cho_Epi, counts = TRUE, data = TRUE, scale.data = FALSE)
seurat2loom(Cho_Epi, "/public/workspace/liuqzh/gastric_cancer/pySCENIC/Fib/Cho_Fib.loom")

####@ parameters
DATABASE="/public/workspace/zhumy/ref/SCENIC"
OUTDIR="/public/workspace/liuqzh/gastric_cancer/pySCENIC/Fib"
# grnboost2
expression_mtx_fname="${OUTDIR}/Cho_Fib.loom"
tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
grnboost2_output="${OUTDIR}/Cho_Fib.adjacencies.tsv"
# ctx
database_fname1=${DATABASE}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
database_fname2=${DATABASE}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
ctx_output="${OUTDIR}/Cho_Fib.regulons.tsv"
# aucell
aucell_output="${OUTDIR}/Cho_Fib.auc_mtx.csv"
aucell_output2="${OUTDIR}/Cho_Fib.auc_mtx.loom"

python ${DATABASE}/arboreto_with_multiprocessing.py \
    $expression_mtx_fname $tfs_fname --method grnboost2 \
    --output $grnboost2_output \
    --num_workers 12 --seed 12345
pyscenic ctx \
    $grnboost2_output \
    $database_fname1 $database_fname2 \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $expression_mtx_fname \
    --output $ctx_output \
    --num_workers 12 --mode "custom_multiprocessing" 
pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 12 --seed 12345
}

##########
{
setwd("/public/workspace/liuqzh/gastric_cancer/pySCENIC/T_cell/")
#########################T cell|NK cell
Cho_Epi=subset(STAD_P,cells=which(STAD_P$celltype=='T cell'|STAD_P$celltype=='NK cell'))
DefaultAssay(Cho_Epi) <- "RNA"

Cho_Epi <- DietSeurat(Cho_Epi, counts = TRUE, data = TRUE, scale.data = FALSE)
seurat2loom(Cho_Epi, "/public/workspace/liuqzh/gastric_cancer/pySCENIC/T_cell/Cho_T_NK.loom")

####@ parameters
DATABASE="/public/workspace/zhumy/ref/SCENIC"
OUTDIR="/public/workspace/liuqzh/gastric_cancer/pySCENIC/T_cell"
# grnboost2
expression_mtx_fname="${OUTDIR}/Cho_T_NK.loom"
tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
grnboost2_output="${OUTDIR}/Cho_T_NK.adjacencies.tsv"
# ctx
database_fname1=${DATABASE}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
database_fname2=${DATABASE}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
ctx_output="${OUTDIR}/Cho_T_NK.regulons.tsv"
# aucell
aucell_output="${OUTDIR}/Cho_T_NK.auc_mtx.csv"
aucell_output2="${OUTDIR}/Cho_T_NK.auc_mtx.loom"

python ${DATABASE}/arboreto_with_multiprocessing.py \
    $expression_mtx_fname $tfs_fname --method grnboost2 \
    --output $grnboost2_output \
    --num_workers 12 --seed 12345
pyscenic ctx \
    $grnboost2_output \
    $database_fname1 $database_fname2 \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $expression_mtx_fname \
    --output $ctx_output \
    --num_workers 12 --mode "custom_multiprocessing" 
pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 12 --seed 12345
}

##########
{
setwd("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Mac/")
#########################Macrophage
Cho_Epi=subset(STAD_P,cells=which(STAD_P$celltype=='Macrophage'))
DefaultAssay(Cho_Epi) <- "RNA"

Cho_Epi <- DietSeurat(Cho_Epi, counts = TRUE, data = TRUE, scale.data = FALSE)
seurat2loom(Cho_Epi, "/public/workspace/liuqzh/gastric_cancer/pySCENIC/Mac/Cho_Mac.loom")

####@ parameters
DATABASE="/public/workspace/zhumy/ref/SCENIC"
OUTDIR="/public/workspace/liuqzh/gastric_cancer/pySCENIC/Mac"
# grnboost2
expression_mtx_fname="${OUTDIR}/Cho_Mac.loom"
tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
grnboost2_output="${OUTDIR}/Cho_Mac.adjacencies.tsv"
# ctx
database_fname1=${DATABASE}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
database_fname2=${DATABASE}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
ctx_output="${OUTDIR}/Cho_Mac.regulons.tsv"
# aucell
aucell_output="${OUTDIR}/Cho_Mac.auc_mtx.csv"
aucell_output2="${OUTDIR}/Cho_Mac.auc_mtx.loom"

python ${DATABASE}/arboreto_with_multiprocessing.py \
    $expression_mtx_fname $tfs_fname --method grnboost2 \
    --output $grnboost2_output \
    --num_workers 12 --seed 12345
pyscenic ctx \
    $grnboost2_output \
    $database_fname1 $database_fname2 \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $expression_mtx_fname \
    --output $ctx_output \
    --num_workers 12 --mode "custom_multiprocessing" 
pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 12 --seed 12345
}

##########
{
setwd("/public/workspace/liuqzh/gastric_cancer/pySCENIC/B_cell/")
#########################Macrophage
Cho_Epi=subset(STAD_P,cells=which(STAD_P$celltype=='B cell'))
DefaultAssay(Cho_Epi) <- "RNA"

Cho_Epi <- DietSeurat(Cho_Epi, counts = TRUE, data = TRUE, scale.data = FALSE)
seurat2loom(Cho_Epi, "/public/workspace/liuqzh/gastric_cancer/pySCENIC/B_cell/Cho_B_cell.loom")

####@ parameters
DATABASE="/public/workspace/zhumy/ref/SCENIC"
OUTDIR="/public/workspace/liuqzh/gastric_cancer/pySCENIC/B_cell"
# grnboost2
expression_mtx_fname="${OUTDIR}/Cho_B_cell.loom"
tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
grnboost2_output="${OUTDIR}/Cho_B_cell.adjacencies.tsv"
# ctx
database_fname1=${DATABASE}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
database_fname2=${DATABASE}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
ctx_output="${OUTDIR}/Cho_B_cell.regulons.tsv"
# aucell
aucell_output="${OUTDIR}/Cho_B_cell.auc_mtx.csv"
aucell_output2="${OUTDIR}/Cho_B_cell.auc_mtx.loom"

python ${DATABASE}/arboreto_with_multiprocessing.py \
    $expression_mtx_fname $tfs_fname --method grnboost2 \
    --output $grnboost2_output \
    --num_workers 12 --seed 12345
pyscenic ctx \
    $grnboost2_output \
    $database_fname1 $database_fname2 \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $expression_mtx_fname \
    --output $ctx_output \
    --num_workers 12 --mode "custom_multiprocessing" 
pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 12 --seed 12345
}

##########Endothelial
{
setwd("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Endothelial/")
#########################Macrophage
Cho_Epi=subset(STAD_P,cells=which(STAD_P$celltype=='Endothelial'))
DefaultAssay(Cho_Epi) <- "RNA"

Cho_Epi <- DietSeurat(Cho_Epi, counts = TRUE, data = TRUE, scale.data = FALSE)
seurat2loom(Cho_Epi, "/public/workspace/liuqzh/gastric_cancer/pySCENIC/Endothelial/Cho_Endothelial.loom")

####@ parameters
DATABASE="/public/workspace/zhumy/ref/SCENIC"
OUTDIR="/public/workspace/liuqzh/gastric_cancer/pySCENIC/Endothelial"
# grnboost2
expression_mtx_fname="${OUTDIR}/Cho_Endothelial.loom"
tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
grnboost2_output="${OUTDIR}/Cho_Endothelial.adjacencies.tsv"
# ctx
database_fname1=${DATABASE}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
database_fname2=${DATABASE}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
ctx_output="${OUTDIR}/Cho_Endothelial.regulons.tsv"
# aucell
aucell_output="${OUTDIR}/Cho_Endothelial.auc_mtx.csv"
aucell_output2="${OUTDIR}/Cho_Endothelial.auc_mtx.loom"

python ${DATABASE}/arboreto_with_multiprocessing.py \
    $expression_mtx_fname $tfs_fname --method grnboost2 \
    --output $grnboost2_output \
    --num_workers 70 --seed 12345
pyscenic ctx \
    $grnboost2_output \
    $database_fname1 $database_fname2 \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $expression_mtx_fname \
    --output $ctx_output \
    --num_workers 70 --mode "custom_multiprocessing" 
pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 70 --seed 12345
}

##########Epithelium
{
setwd("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/")
#########################Macrophage
Cho_Epi=subset(STAD_P,cells=which(STAD_P$celltype=='Epithelium'))
DefaultAssay(Cho_Epi) <- "RNA"

Cho_Epi <- DietSeurat(Cho_Epi, counts = TRUE, data = TRUE, scale.data = FALSE)
seurat2loom(Cho_Epi, "/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium/Cho_Epithelium.loom")

####@ parameters
DATABASE="/public/workspace/zhumy/ref/SCENIC"
OUTDIR="/public/workspace/liuqzh/gastric_cancer/pySCENIC/Epithelium"
# grnboost2
expression_mtx_fname="${OUTDIR}/Cho_Epithelium.loom"
tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
grnboost2_output="${OUTDIR}/Cho_Epithelium.adjacencies.tsv"
# ctx
database_fname1=${DATABASE}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
database_fname2=${DATABASE}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
ctx_output="${OUTDIR}/Cho_Epithelium.regulons.tsv"
# aucell
aucell_output="${OUTDIR}/Cho_Epithelium.auc_mtx.csv"
aucell_output2="${OUTDIR}/Cho_Epithelium.auc_mtx.loom"

python ${DATABASE}/arboreto_with_multiprocessing.py \
    $expression_mtx_fname $tfs_fname --method grnboost2 \
    --output $grnboost2_output \
    --num_workers 90 --seed 12345
pyscenic ctx \
    $grnboost2_output \
    $database_fname1 $database_fname2 \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $expression_mtx_fname \
    --output $ctx_output \
    --num_workers 90 --mode "custom_multiprocessing" 
pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 90 --seed 12345
}

##########Mast cell
{
setwd("/public/workspace/liuqzh/gastric_cancer/pySCENIC/Mast_cell/")
#########################Macrophage
Cho_Epi=subset(STAD_P,cells=which(STAD_P$celltype=='Mast cell'))
DefaultAssay(Cho_Epi) <- "RNA"

Cho_Epi <- DietSeurat(Cho_Epi, counts = TRUE, data = TRUE, scale.data = FALSE)
seurat2loom(Cho_Epi, "/public/workspace/liuqzh/gastric_cancer/pySCENIC/Mast_cell/Cho_Mast_cell.loom")

####@ parameters
DATABASE="/public/workspace/zhumy/ref/SCENIC"
OUTDIR="/public/workspace/liuqzh/gastric_cancer/pySCENIC/Mast_cell"
# grnboost2
expression_mtx_fname="${OUTDIR}/Cho_Mast_cell.loom"
tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
grnboost2_output="${OUTDIR}/Cho_Mast_cell.adjacencies.tsv"
# ctx
database_fname1=${DATABASE}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
database_fname2=${DATABASE}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
ctx_output="${OUTDIR}/Cho_Epithelium.regulons.tsv"
# aucell
aucell_output="${OUTDIR}/Cho_Mast_cell.auc_mtx.csv"
aucell_output2="${OUTDIR}/Cho_Mast_cell.auc_mtx.loom"

python ${DATABASE}/arboreto_with_multiprocessing.py \
    $expression_mtx_fname $tfs_fname --method grnboost2 \
    --output $grnboost2_output \
    --num_workers 90 --seed 12345
pyscenic ctx \
    $grnboost2_output \
    $database_fname1 $database_fname2 \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $expression_mtx_fname \
    --output $ctx_output \
    --num_workers 90 --mode "custom_multiprocessing" 
pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 90 --seed 12345
}









