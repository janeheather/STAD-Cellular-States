###########################
bytlib load R-4.0.2
bytlib load gcc
R

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)

All_MP_subtype_Result=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv',header=T,sep=',')

# 假设每个数据框的第一列是一个可以用于对齐数据的唯一标识符
# 使用这个标识符作为合并的基础
Com2=c('B_MP_5','T_cell_MP_1','Mast_MP_3','Endo_MP_6','Maligant_MP_1','T_cell_MP_3','B_MP_4','T_cell_MP_4','B_MP_3','Mac_MP_2','Endo_MP_3')

merged_data=All_MP_subtype_Result[,Com2]

library(igraph)

data=merged_data
# 创建节点数据框
nodes <- data.frame(name = colnames(data))
nodes$size <- sapply(data, function(x) sum(!is.na(x)))  # 计算每列中非空元素的数量

# 创建边数据框
edges <- expand.grid(from = colnames(data), to = colnames(data))
edges <- subset(edges, from != to)  # 移除自环

# 计算杰卡德相似性
jaccard_similarity <- function(col1, col2, data) {
    set1 <- data[[col1]][!is.na(data[[col1]])]
    set2 <- data[[col2]][!is.na(data[[col2]])]
    return(length(intersect(set1, set2)) / length(union(set1, set2)))
}

edges$weight <- apply(edges, 1, function(x) jaccard_similarity(x['from'], x['to'], data))

# 筛选权重大于0.02的边
edges <- subset(edges, weight > 0.02)

# 创建图对象
g <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)

# 简化图对象以删除任何重复的边或自环
g <- simplify(g)

# 设置节点颜色
V(g)$color <- ifelse(V(g)$name == "Maligant_MP_1", "#D1832C",
                     ifelse(V(g)$name == "EMT", "#5B93D4",
                            ifelse(V(g)$name == "NLT", "#A97598", "skyblue")))
V(g)$color <- ifelse(grepl("Fib", V(g)$name), "#509285", V(g)$color)
V(g)$color <- ifelse(grepl("Mac", V(g)$name), "#701750", V(g)$color)
V(g)$color <- ifelse(grepl("T_cell", V(g)$name), "#9F70A7", V(g)$color)
V(g)$color <- ifelse(grepl("B", V(g)$name), "#A2855B", V(g)$color)
V(g)$color <- ifelse(grepl("Mast_cell", V(g)$name), "#427E39", V(g)$color)
V(g)$color <- ifelse(grepl("Endo", V(g)$name), "#5891D2", V(g)$color)

# 设置节点大小
V(g)$size <- nodes$size / 10  # 调整大小系数以便于观察

# 绘图
p = plot(g, 
		vertex.size=V(g)$size, 
		vertex.color=V(g)$color, 
		edge.width=E(g)$weight*30,  # 调整以便于观察
		edge.color="grey",
		main="Network Graph",
		layout=layout_with_fr(g))

pdf('/public/workspace/liuqzh/gastric_cancer/Com/Jaccard_network_plot_Comaa.pdf',height=6,width=6)
plot(g, 
		vertex.size=V(g)$size, 
		vertex.color=V(g)$color, 
		edge.width=E(g)$weight*10,  # 调整以便于观察
		edge.color="grey",
		main="Network Graph",
		layout=layout_with_fr(g))
dev.off()

#############
#############
###########################
bytlib load R-4.0.2
bytlib load gcc
R

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)

All_MP_subtype_Result=read.csv('/public/workspace/liuqzh/gastric_cancer/scalop-master/All_MP_subtype_Result.csv',header=T,sep=',')

# 假设每个数据框的第一列是一个可以用于对齐数据的唯一标识符
# 使用这个标识符作为合并的基础
Com4=c('Endo_MP_2','Mast_MP_2','Fib_MP_3','T_cell_MP_5','Maligant_MP_3','Mast_MP_1','T_cell_MP_6','B_MP_2','T_cell_MP_2','B_MP_1','Endo_MP_4')

merged_data=All_MP_subtype_Result[,Com4]

library(igraph)

data=merged_data
# 创建节点数据框
nodes <- data.frame(name = colnames(data))
nodes$size <- sapply(data, function(x) sum(!is.na(x)))  # 计算每列中非空元素的数量

# 创建边数据框
edges <- expand.grid(from = colnames(data), to = colnames(data))
edges <- subset(edges, from != to)  # 移除自环

# 计算杰卡德相似性
jaccard_similarity <- function(col1, col2, data) {
    set1 <- data[[col1]][!is.na(data[[col1]])]
    set2 <- data[[col2]][!is.na(data[[col2]])]
    return(length(intersect(set1, set2)) / length(union(set1, set2)))
}

edges$weight <- apply(edges, 1, function(x) jaccard_similarity(x['from'], x['to'], data))

# 筛选权重大于0.02的边
edges <- subset(edges, weight > 0.02)

# 创建图对象
g <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)

# 简化图对象以删除任何重复的边或自环
g <- simplify(g)

# 设置节点颜色
V(g)$color <- ifelse(V(g)$name == "Maligant_MP_1", "#D1832C",
                     ifelse(V(g)$name == "EMT", "#5B93D4",
                            ifelse(V(g)$name == "Maligant_MP_3", "#A97598", "skyblue")))
V(g)$color <- ifelse(grepl("Fib", V(g)$name), "#509285", V(g)$color)
V(g)$color <- ifelse(grepl("Mac", V(g)$name), "#701750", V(g)$color)
V(g)$color <- ifelse(grepl("T_cell", V(g)$name), "#9F70A7", V(g)$color)
V(g)$color <- ifelse(grepl("B", V(g)$name), "#A2855B", V(g)$color)
V(g)$color <- ifelse(grepl("Mast_cell", V(g)$name), "#427E39", V(g)$color)
V(g)$color <- ifelse(grepl("Endo", V(g)$name), "#5891D2", V(g)$color)

# 设置节点大小
V(g)$size <- nodes$size / 10  # 调整大小系数以便于观察

# 绘图
p = plot(g, 
		vertex.size=V(g)$size, 
		vertex.color=V(g)$color, 
		edge.width=E(g)$weight*30,  # 调整以便于观察
		edge.color="grey",
		main="Network Graph",
		layout=layout_with_fr(g))

pdf('/public/workspace/liuqzh/gastric_cancer/Com/Jaccard_network_plot_Com4.pdf',height=6,width=6)
plot(g, 
		vertex.size=V(g)$size, 
		vertex.color=V(g)$color, 
		edge.width=E(g)$weight*10,  # 调整以便于观察
		edge.color="grey",
		main="Network Graph",
		layout=layout_with_fr(g))
dev.off()