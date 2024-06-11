########
bytlib load R-4.0.2
bytlib load gcc
R
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(estimate)

Maligant_MP_1<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_1.txt',header=T,row.names=1)
Maligant_MP_2<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_2.txt',header=T,row.names=1)
Maligant_MP_3<-read.table('/public/workspace/liuqzh/gastric_cancer/Simplicity/Maligant_MP_3.txt',header=T,row.names=1)

normalize <- function(x) {
return((x - min(x)) / (max(x) - min(x)))
}
Maligant_MP_1=arrange(Maligant_MP_1,desc(Simplicity))
Maligant_MP_1$Simplicity_score=normalize(Maligant_MP_1$Simplicity)
Maligant_MP_2=arrange(Maligant_MP_2,desc(Simplicity))
Maligant_MP_2$Simplicity_score=normalize(Maligant_MP_2$Simplicity)
Maligant_MP_3=arrange(Maligant_MP_3,desc(Simplicity))
Maligant_MP_3$Simplicity_score=normalize(Maligant_MP_3$Simplicity)

colnames(Maligant_MP_1)[15]='ADDS'
colnames(Maligant_MP_2)[15]='ADDS'
colnames(Maligant_MP_3)[15]='ADDS'

DATA_MP=rbind(Maligant_MP_1,Maligant_MP_2,Maligant_MP_3)
DATA_MP$MP1 = -log10(DATA_MP$Maligant_MP_1_pval)
DATA_MP$MP2 = -log10(DATA_MP$Maligant_MP_2_pval)
DATA_MP$MP3 = -log10(DATA_MP$Maligant_MP_3_pval)

DATA_MP$simplicity='Low'
DATA_MP$simplicity[which(DATA_MP$Simplicity_score >= 0.8)]='High'

tmb=read.table("/public/workspace/liuqzh/gastric_cancer/TCGA_TMB/tmb.txt",header=T,row.names=1)
rownames(tmb) <- substr(rownames(tmb), 1, 16)

tmb$sample=rownames(tmb)
DATA_MP$sample=rownames(DATA_MP)
DATA_MP_TMB=merge(DATA_MP,tmb,'sample')

write.table(DATA_MP_TMB, "/public/workspace/liuqzh/gastric_cancer/TCGA_TMB/DATA_MP_TMB.txt", row.names = F, sep = "\t")

library(ggplot2)
library(ggpubr)

# 假设您的数据框叫 df，它包含的列有 Group、total_perMB_log 和 total_perMB
# 确保替换下面的 df 为您的实际数据框变量名
df=DATA_MP_TMB

df$Subtype <- factor(df$Subtype, levels = c("MP_1", "MP_2", "MP_3"))  # 设置分组的顺序

# 生成箱形图并添加散点图层
p <- ggplot(df, aes(x = Subtype, y = total_perMB_log, color = Subtype)) +
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
  labs(title = "TMB with STAD samples by Subtype", x = "Subtype", y = "total_perMB_log")  # 添加标签

# 计算并添加指定组间的 p 值
p <- p + stat_compare_means(comparisons = list(c("MP_1", "MP_2"),c("MP_2", "MP_3"),c("MP_1", "MP_3")), 
                            label = "p.signif", label.y = c(2.5, 3, 3.5))

pdf('/public/workspace/liuqzh/gastric_cancer/Simplicity/TCGA_TMB.pdf',height=6,width=4)
print(p)
dev.off()

# 打印图形
print(p)