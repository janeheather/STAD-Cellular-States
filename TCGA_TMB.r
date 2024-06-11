
bytlib load R-4.0.2
bytlib load gcc
R

library(dplyr)
library(maftools)
all.maf <- list.files(path = "/public/workspace/yumiao/STAD/TCGA-STAD/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/", pattern = ".gz",full.names = T, recursive = T)#full.names = T返回绝对路径
# 看看前3个
all.maf[1:3]

# 读取所有文件
maf.list <- lapply(all.maf, data.table::fread,#read.table, 
                   skip = 7, 
                   sep = "\t", 
                   header = T)
# 查看列数是不是一样 140列
lapply(maf.list, dim)

# 合并
maf.merge <- do.call(rbind,maf.list)#194729行

#读取maf文件
library(maftools)
maf1 <- read.maf(maf.merge)
tmb <- tmb(maf = maf1,
           captureSize = 50,
           logScale = T)
head(tmb)

tmb=as.data.frame(tmb)

write.table(tmb, "/public/workspace/liuqzh/gastric_cancer/TCGA_TMB/tmb.txt", row.names = F, sep = "\t")

#########
tmb=read.table("/public/workspace/liuqzh/gastric_cancer/TCGA_TMB/tmb.txt",header=T,row.names=1)
Subtype_III_tcga=read.csv("/public/workspace/liuqzh/gastric_cancer/cibersort/CIBERSORTx_Job9_Results.txt",header=T,sep='\t',row.names=1)

rownames(tmb) <- substr(rownames(tmb), 1, 12)

tmb_tmp=tmb[intersect(rownames(Subtype_III_tcga),rownames(tmb)),]
Subtype_III_tcga=Subtype_III_tcga[intersect(rownames(tmb_tmp),rownames(Subtype_III_tcga)),]

Subtype_III_tcga_TMB=cbind(Subtype_III_tcga,tmb_tmp)

write.table(Subtype_III_tcga_TMB, "/public/workspace/liuqzh/gastric_cancer/cibersort/Subtype_III_tcga_TMB.txt", row.names = F, sep = "\t")

Subtype_III_tcga_TMB$Subtype='MAX'
Subtype_III_tcga_TMB$Subtype[which(Subtype_III_tcga_TMB$MP1 > 0.5)]='MP1'
Subtype_III_tcga_TMB$Subtype[which(Subtype_III_tcga_TMB$MP2 > 0.5)]='MP2+MP3'
Subtype_III_tcga_TMB$Subtype[which(Subtype_III_tcga_TMB$MP3 > 0.5)]='MP2+MP3'

library(ggplot2)
library(ggpubr)

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
df=Subtype_III_tcga_TMB
length(which(Subtype_III_tcga_TMB$Subtype == 'MAX'))
df=df[which(df$Subtype != 'MAX'),]

df$group <- factor(df$Subtype, levels = c("MP1", "MP2+MP3"))  # 设置分组的顺序
df$Color <- ifelse(df$total_perMB > 10, "red", "lightblue")  # 设置颜色条件

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = group, y = total_perMB_log, color = Color)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +  # 不显示离群值，且箱形图无填充
  geom_jitter(size = 2, width = 0.2) +  # 添加散点图层
  scale_color_manual(values = c("red" = "red", "lightblue" = "lightblue")) +  # 设置颜色
  theme_minimal() +  # 使用简洁主题
  theme(  # 调整主题去掉不必要的元素
    panel.grid.major = element_blank(),  # 移除网格线
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),  # 移除背景
    axis.line = element_line(color = "black"),  # 坐标轴线
    legend.position = "none"  # 不显示图例
  ) +
  labs(title = "TMB with STAD samples by Group", x = "Group", y = "total_perMB_log")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP1", "MP2+MP3")), 
                            label = "p.signif", label.y = c(2.2, 2.5, 2.8))

pdf('/public/workspace/liuqzh/gastric_cancer/TCGA_TMB/TCGA_TMB_MP1_MP2_MP3.pdf',height=6,width=4)
print(p)
dev.off()

# 打印图形
print(p)
#########

head(df)
data_f = df[which(df$MP1 > 0.5),]

p1=ggplot(data=data_f, aes(x=MP1, y=total_perMB_log))+geom_point(color="red")+stat_smooth(method="lm",se=TRUE)+stat_cor(data=data_f, method = "pearson")
# 你的数据框名为data_f，确保已正确加载
p1<-ggplot(data = data_f, aes(x = MP1, y = total_perMB_log)) +
	geom_point(color = "red") +  # 添加红色的点
	stat_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加线性模型拟合线，并显示置信区间，颜色设为蓝色
	stat_cor(method = "pearson", label.x = 0.5, label.y = 2) +  # 添加皮尔逊相关系数
	theme_minimal() +  # 使用简洁主题
	labs(
		title = "Relationship between MP1 and total_perMB_log",
		x = "MP1",
		y = "Total per MB log"
	) +  # 添加图形标题和轴标题
	theme(
		plot.title = element_text(hjust = 0.5),  # 标题居中
		axis.text = element_text(color = "gray20"),  # 轴文字颜色
		axis.title = element_text(color = "gray20"),  # 轴标题颜色
		panel.grid.major = element_blank(),  # 移除主要网格线
		panel.grid.minor = element_blank(),  # 移除次要网格线
		panel.background = element_rect(fill = "white", color = "gray50")  # 背景颜色
	)


pdf('/public/workspace/liuqzh/gastric_cancer/TCGA_TMB/TCGA_TMB_MP1_corr.pdf',height=6,width=6)
print(p1)
dev.off()

ggplot(data=df, aes(x=MP2, y=total_perMB_log))+geom_point(color="red")

cor.test(df$MP1,df$total_perMB_log)

