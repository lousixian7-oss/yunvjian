##免费的午餐limma+火山图代码##


#设置环境变量
library(tidyverse)
library(limma)
exp<-read.csv("去除批次效应后表达矩阵.csv",row.names = 1)#读入表达矩阵文件

##数据清洗###
summary(as.numeric(exp))

# 均数不一致需标准化 ------------------------------
boxplot(exp,outline=FALSE, notch=F , las=2)


##去除低表达的基因 ------------------------------

exp = exp[rowMeans(exp)>1,] #根据自己的需要去除低表达基因，也可以卡其它阈值

#读入分类文件------------------------------

fen<-read.csv("分类（新）.csv",row.names = 1)#读入表达矩阵表头的分类文件

##开始进行差异分析-----------------------------

group_list <- factor(fen$group,levels = c("NC","PD"))
design <- model.matrix(~group_list)
con <- lmFit(exp,design)
con2 <- eBayes(con)
DEG1 <- topTable(con2, coef = 2, number = Inf)
DEG2 = na.omit(DEG1) 

##导出差异分析结果-----------------------------

write.csv(DEG2,"差异分析结果11.csv")


#筛选差异基因----------------------------

logFC_cut = 1
p_cut = 0.05

type1 = (DEG2$P.Value < p_cut)&(DEG2$logFC < -logFC_cut)
type2 = (DEG2$P.Value < p_cut)&(DEG2$logFC > logFC_cut)

DEG2$type = ifelse(type1,"Down",ifelse(type2,"Up","NOT"))

head(DEG2)

write.csv(DEG2,"差异分析标注11.csv")


##绘制火山图-------------------带标记基因---------

library(ggplot2)
library(ggrepel)
library(dplyr)


#===============================
# DEG分类
#===============================

logFC_cut <- 0.25
p_cut <- 0.05


DEG2$type <- ifelse(
  DEG2$adj.P.Val < p_cut &
    DEG2$logFC > logFC_cut,
  "Up",
  ifelse(
    DEG2$adj.P.Val < p_cut &
      DEG2$logFC < -logFC_cut,
    "Down",
    "NS"
  )
)


#===============================
# 指定标记基因BRCA1
# AGT
# C3
# CYP3A4
# CXCR4
# FOS

#===============================

label_genes <- c(
  "BRCA1",
  "AGT",
  "C3",
  "FOS",
  "CYP3A4","CXCR4","PI16","COL1A1","FN1"
)


# 把行名转成基因列
DEG2$gene <- rownames(DEG2)



label_df <- DEG2 %>%
  filter(gene %in% label_genes)



#===============================
# SCI volcano
#===============================


p <- ggplot(
  DEG2,
  aes(
    x = logFC,
    y = -log10(adj.P.Val)
  )
)+
  
  geom_point(
    aes(color=type),
    size=1.5,
    alpha=0.7
  )+
  
  scale_color_manual(
    values=c(
      "Down"="#00BFC4",
      "NS"="grey80",
      "Up"="#F8766D"
    )
  )+
  
  geom_hline(
    yintercept=-log10(0.05),
    linetype="dashed",
    color="grey40"
  )+
  
  geom_vline(
    xintercept=c(-1,1),
    linetype="dashed",
    color="grey40"
  )+
  
  geom_label_repel(
    data=label_df,
    aes(label=gene),
    size=4,
    box.padding=0.5,
    max.overlaps=20,
    color="black",
    fill="white"
  )+
  
  theme_classic(base_size=14)+
  
  labs(
    x="log2 Fold Change",
    y="-log10(FDR)"
  )+
  
  theme(
    legend.title=element_blank(),
    legend.position="right",
    axis.title=element_text(face="bold"),
    axis.text=element_text(color="black")
  )


p

library(pheatmap)
library(RColorBrewer)

# DEG排序
DEG_sig <- DEG2 %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val)

top_gene <- rownames(DEG_sig)[1:50]


heat_exp <- exp[top_gene, ]

# Z-score
heat_exp <- t(scale(t(heat_exp)))


# annotation
annotation_col <- data.frame(
  Group = fen$group,
  Dataset = fen$dataset
)

rownames(annotation_col) <- colnames(exp)


ann_colors <- list(
  Group=c(
    NC="#00BFC4",
    PD="#F8766D"
  ),
  Dataset=c(
    GSE10334="#E7B800",
    GSE156993="#FC8D62",
    GSE23586="#96CEB4"
  )
)


pheatmap(
  heat_exp,
  annotation_col=annotation_col,
  annotation_colors=ann_colors,
  cluster_cols=TRUE,
  cluster_rows=TRUE,
  show_colnames=FALSE,
  fontsize_row=8,
  color=colorRampPalette(
    c("#3B4CC0","white","#B40426")
  )(100)
)

ggsave(
  "Bulk_DEG_volcano_TXNRD1.pdf",
  p,
  width=6,
  height=5
)