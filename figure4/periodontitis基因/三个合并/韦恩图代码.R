# 加载必要的包
library(ggvenn)
library(eulerr)
library(scales)
library(dplyr)  

# 读取数据
data <- read.csv("pd并集基因.csv")

clean_genes <- function(x) {
  x <- as.character(x)
  x <- trimws(x)                     # 去除首尾空格
  x <- x[!is.na(x) & x != ""]        # 去掉 NA 和空字符串
  return(unique(x))                  # 去重
}
# 提取两列基因作为集合
comgene <- list(
  'Periodontitis' = clean_genes(data$X1),
  'Yunvjian' = clean_genes(data$X2),
  'DEGs' = clean_genes(data$X3)
)


# 计算交集基因
com_genes <- intersect(comgene$Periodontitis, comgene$Yunvjian)
com_genes <- intersect(com_genes, comgene$DEGs)

com_genes<- clean_genes(com_genes)
#保存为 CSV
{
write.csv(data.frame(Intersection_Genes = com_genes),
          "pd&yunajian交集基因.csv", row.names = FALSE)
}

#普通韦恩图
p1 <- ggvenn(
  comgene,
  show_percentage = F,
  show_elements = FALSE,
  label_sep = ",",
  digits = 2,
  stroke_color = NA,  # 去除边框
  fill_color = c("#6699CC", "#9966CC","#669966"), 
  set_name_color = c("#6699CC", "#9966CC","#669966"),
  text_color = "black",  # 调整文字颜色
  text_size = 5          # 增大文字
)
print(p1)

# 绘制比例图
p2 <- euler(comgene)
plot(
  p2,
  labels = list(col = "black", font = 1, cex = 1.5),  # 调整字体大小和颜色
  edges = NULL,  # 去除边框
  fills = c("#C2E0F7","#FFD39B"),  
  alpha = 0.5,
  quantities = list(cex = 1.5, col = 'black')  # 设置数量文字更突出
)

