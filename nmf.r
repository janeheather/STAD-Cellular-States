bytlib load R-4.0.2
bytlib load gcc
R

library(reshape2)
library(NMF)
library(ggplot2)
library(scales)

# Custom color palette
library(RColorBrewer)
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

# 用于轮廓系数的计算
library(cluster)    
library(dplyr)

# -------------------------------------------------------------------------------------------
# Function for selecting robust nonnegative matrix factorization (NMF) programs
# ------------------------------------------------------------------------------------------- 

robust_nmf_programs <- function(nmf_programs, intra_min = 60, intra_max = 30, inter_filter=T, inter_min = 20) {
	# Select NMF programs based on the minimum overlap with other NMF programs from the sample
	intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y)))))
	intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))
	
	###模块样本稳定性，在样本内部和其他模块至少要有75的相同基因，确定之后固定Module
	nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min])
	names(nmf_sel) <- names(nmf_programs)

	# Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
	# ii) the minimum overlap with programs from another cell line
	
	###每个Module间的相似性
	nmf_sel_unlist <- do.call(cbind, nmf_sel)
	inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs

	
	final_filter <- NULL 
	for(i in names(nmf_sel)) {
		###计算某个样本内部的Module与其他样本的全部Module的相似性
		a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
		b <- sort(apply(a, 2, max), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
		
		###计算某个样本内部的Module与其他样本的全部Module的相似性，选出最高的相似性，如果大于30个基因Module保留
		###样本的特异性，特异性太高的Module删除
		if(inter_filter==T) b <- b[b>=inter_min] # selects programs with a maximum intersection of at least 40
			if(length(b) > 1) {
				c <- names(b[1]) 
				for(y in 2:length(b)) {
					###计算样本剩下Module之间的相似性，如果相似性太高去掉，小于20个基因
					if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 20 with a previously selected programs
				}
				final_filter <- c(final_filter, c)
			} 
			else {
				final_filter <- c(final_filter, names(b))
		}
	}
	return(final_filter)                                                      
}

nmf_score_result=readRDS('/public/workspace/liuqzh/gastric_cancer/NMF_DATA/nmf_score_result_0.4_0.7_0.2_DATA.rds')
nmf_programs=readRDS('/public/workspace/liuqzh/gastric_cancer/NMF_DATA/nmf_programs_0.4_0.7_0.2_DATA.rds')
##

Genes_nmf_w_basis=nmf_score_result

## Parameters 
intra_min_parameter <- 70
intra_max_parameter <- 30
inter_min_parameter <- 20

# for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
nmf_filter       <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter=T, inter_min = inter_min_parameter)  
mmf=nmf_filter

nmf_programs     <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter),drop=F])
nmf_programs     <- do.call(cbind, nmf_programs)

# calculate similarity between programs
nmf_intersect         <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc     <- hclust(dist(100-nmf_intersect,method='euclidean'), method="complete")

load(file = '/public/workspace/liuqzh/gastric_cancer/NMF_DATA/NMF_75_30_20.RData')
# 计算不同数量的聚类对应的轮廓系数
max_k <- 10  # 可以根据需要测试的最大聚类数来调整
sil_width <- numeric(max_k - 1)
for (k in 2:max_k) {
	# 切割树状图以获得k个聚类
	cluster_labels <- cutree(nmf_intersect_hc, k)
	# 计算轮廓系数
	sil_width[k - 1] <- silhouette(cluster_labels, as.dist(100 - nmf_intersect))[, 3] %>% mean()
}

# 创建一个数据框来存储聚类数和对应的轮廓宽度
sil_df <- data.frame(k = 2:max_k, sil_width)

# 绘制轮廓系数图
pdf('/public/workspace/liuqzh/gastric_cancer/NMF_DATA/Silhouette.pdf',width=5,height=5,useDingbats=F)
ggplot(sil_df, aes(x = k, y = sil_width)) +
		geom_line() +
		geom_point() +
		scale_x_continuous(breaks = 2:max_k) +
		labs(x = "Number of clusters", y = "Average silhouette width") +
		ggtitle("Silhouette Coefficient")+
		theme_minimal()+
	theme(
		plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 调整标题样式
		axis.title = element_text(size = 14, face = "bold"),  # 调整轴标题样式
		axis.text = element_text(size = 12),  # 调整轴文本大小
		legend.position = "none"  # 移除图例（如果不必要）
	) +
		scale_colour_manual(values = c("#00BFC4", "#F8766D"))  # 自定义颜色（根据需要调整）
dev.off()

k_best=which(sil_df$sil_width %in% max(sil_df$sil_width))+1

# 获得轮廓系数为时的聚类结果
cluster_labels <- cutree(nmf_intersect_hc, k=k_best)

pdf('/public/workspace/liuqzh/mf_intersect_hc.pdf',width=50,height=10,useDingbats=F)
plot(nmf_intersect_hc)
dev.off()

cluster_labels=as.matrix(cluster_labels)
cluster_labels=cluster_labels[intersect(rownames(nmf_intersect),rownames(cluster_labels)),]

#save.image(file = '/public/workspace/liuqzh/gastric_cancer/NMF_DATA/NMF_75_30_20.RData')
### Parameters for clustering

bytlib load R-4.0.2
bytlib load gcc
R

library(reshape2)
library(NMF)
library(ggplot2)
library(scales)

# Custom color palette
library(RColorBrewer)
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

# 用于轮廓系数的计算
library(cluster)    
library(dplyr)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
load(file = '/public/workspace/liuqzh/gastric_cancer/NMF_DATA/NMF_75_30_20.RData')
#write.table(nmf_intersect,'/public/workspace/liuqzh/gastric_cancer/NMF_DATA/nmf_intersect_data.txt',col.names=T, sep='\t')
#write.table(cluster_labels,'/public/workspace/liuqzh/gastric_cancer/NMF_DATA/cluster_labels_data.txt',col.names=T, sep='\t')

pdf('/public/workspace/liuqzh/gastric_cancer/NMF_DATA/Cluster_Dendrogram.pdf',width=60,height=8,useDingbats=F)
plot(nmf_intersect_hc)
dev.off()

cluster_labels=as.data.frame(cluster_labels)
# 将聚类标签转换为因子，以便在注释中使用
colnames(cluster_labels)='cluster_labels'
cluster_labels$cluster_labels <- as.factor(cluster_labels$cluster_labels)

# 构建一个注释数据框，用于热图的列注释
annotation_col <- data.frame(Cluster = cluster_labels$cluster_labels)
rownames(annotation_col) <- rownames(cluster_labels)

# 你的自定义 magma 配色方案
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
breaks <- seq(0, 20, length.out = length(custom_magma))

pdf('/public/workspace/liuqzh/gastric_cancer/NMF_DATA/Heatmap_subtype.pdf',width=8,height=8,useDingbats=F)
pheatmap(nmf_intersect,
         color = custom_magma,  # 使用你的自定义配色方案
         breaks = breaks,  # 设置颜色的断点
         annotation_col = annotation_col,
		 cluster_rows = TRUE,  # 假设你不想对行进行聚类
         cluster_cols = TRUE,  # 假设你不想对列进行聚类
         show_rownames = FALSE,  # 显示行名
         show_colnames = FALSE,  # 显示列名
         border_color = NA,  # 移除单元格边界
         legend = TRUE  # 显示图例
)
dev.off()


nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]


nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect) 

ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
	geom_tile() + 
	scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
	scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
	theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
	theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
	theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
	guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

# ----------------------------------------------------------------------------------------------------
# Cluster selected NMF programs to generate MPs
# ----------------------------------------------------------------------------------------------------

save.image(file = '/public/workspace/liuqzh/gastric_cancer/NMF_DATA/NMF_75_30_20.RData')

load(file = '/public/workspace/liuqzh/gastric_cancer/NMF_DATA/NMF_75_30_20.RData')

MP_list <- c()

for(j in 1:3){
	nmf_programs_tmp=nmf_programs[,intersect(colnames(nmf_programs),names(cluster_labels[which(cluster_labels==j)]))]
	nmf_tmp <- apply(nmf_programs_tmp,2,function(x) apply(nmf_programs_tmp,2,function(y) length(intersect(x,y))))

	Sorted_intersection = as.matrix(sort(apply(nmf_tmp,2,function(x) (length(which(x>=ceiling(0.2*ncol(nmf_programs_tmp))))-1)),decreasing = TRUE))
	Sorted_intersection = Sorted_intersection[which(Sorted_intersection > 0),]
	
	nmf_programs_tmp = nmf_programs_tmp[,which(colnames(nmf_programs_tmp) %in% names(Sorted_intersection))]
	
	Genes_MP = names(sort(table(nmf_programs_tmp),decreasing = TRUE)[1:100])
	
	MP_list[[paste0("MP_",j)]] = Genes_MP
}
MP_list <-  do.call(cbind, MP_list)

MP_list=as.data.frame(MP_list)
MP_list$MP_1=gsub('\\.','-',MP_list$MP_1)
MP_list$MP_2=gsub('\\.','-',MP_list$MP_2)
MP_list$MP_3=gsub('\\.','-',MP_list$MP_3)


MP_list$MP_4=gsub('\\.','-',MP_list$MP_4)
MP_list$MP_5=gsub('\\.','-',MP_list$MP_5)
MP_list$MP_6=gsub('\\.','-',MP_list$MP_6)	