# 加载必要的包
library(ggvenn)
library(eulerr)
library(scales)
library(dplyr)  

# 读取数据
data <- read.csv("输入数据.csv")

# 提取两列基因作为集合
comgene <- list(
  'Disgenet' = unique(data$Group1),
  'Drugbank' = unique(data$Group2),
  'Genecard' = unique(data$Group3),
  'eQTL' = unique(data$eqtl),
  'WGCNA' = unique(data$GEO)
)


# 计算交集基因
com_genes <- union(comgene$Disgenet, comgene$Drugbank)
com_genes <- union(com_genes,comgene$Genecard)
com_genes <- union(com_genes,comgene$eQTL)
com_genes <- union(com_genes,comgene$WGCNA)
com_genes <-unique(com_genes)
#保存为 CSV
{
write.csv(data.frame(Intersection_Genes = com_genes),
          "并集基因.csv", row.names = FALSE)
}
library(venn)  
p1 <-venn(comgene)
#普通韦恩图
p1 <- ggvenn(
  comgene,
  show_percentage = F,
  show_elements = FALSE,
  label_sep = ",",
  digits = 2,
  stroke_color = NA,  # 去除边框
  fill_color = c("#C2E0F7", "#A4CBA8","#FFD39B","#FF9999","#CC99FF"), 
  set_name_color = c("#C2E0F7", "#A4CBA8","#FFD39B","#FF9999","#CC99FF"),
  text_color = "black",  # 调整文字颜色
  text_size = 5          # 增大文字
)
print(p1)

# 绘制比例图
p2 <- euler(comgene)
plot(
  p2,
  labels = list(col = "black", font = 1, cex = 1.0),  # 调整字体大小和颜色
  edges = NULL,  # 去除边框
  fills = c("#C2E0F7", "#A4CBA8","#FFD39B","#FF9999","#CC99FF"),  
  alpha = 0.5,
  quantities = list(cex = 1.0, col = 'black')  # 设置数量文字更突出
)

