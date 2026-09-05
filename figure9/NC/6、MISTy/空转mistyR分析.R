library(Seurat)
library(dplyr)
library(tibble)
library(mistyR)
library(ggplot2)
library(patchwork)
library(viridis)
library(tidyr)

set.seed(123)


##  读取数据
load("牙龈NC_注释.rda")  
## 查看空间切片，确认空间对象是否正常
SpatialPlot(stRNA)

## 读取 RCTD 反卷积结果
load("myRCTD合并.rda")  
myRCTD<-readRDS("stRNA_RCTD_refinedFib_ECM.rds")
##  提取 RCTD 细胞比例矩阵
rctd_mat <- as.data.frame(myRCTD@results$weights)

{## 去掉全 0 的 spot，如果某个 spot 所有细胞类型比例都是 0，说明该 spot 没有有效反卷积结果，应该删除
  rctd_mat <- rctd_mat[rowSums(rctd_mat) > 0, ]
  
  ## 将 RCTD 权重转换为比例
  rctd_mat <- rctd_mat / rowSums(rctd_mat)
  
  ## 将 NA 替换为 0，NA 可能来自除以 0 或缺失值，前面已经去掉全 0 spot，这里只是保险处理
  rctd_mat[is.na(rctd_mat)] <- 0
  
  ## 去掉没有变化的细胞类型
  rctd_mat <- rctd_mat[, apply(rctd_mat, 2, sd) > 0]
}


##  提取空间坐标
img <- Images(stRNA)[1]#Images(stRNA)[1] 获取第一个空间切片图像名称

## 从 Seurat 空间对象中提取 tissue coordinates
## 行名通常是 spot barcode
coords <- GetTissueCoordinates(stRNA, image = img)

coords <- as.data.frame(coords)

{
  ## 不同 Seurat 版本中，空间坐标列名可能不同：1）imagecol / imagerow，2）x / y， 3）col / row
  ## 下面代码用于自动识别并统一命名为 x 和 y
  if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    
    ## Seurat 常见坐标列
    coords <- coords[, c("imagecol", "imagerow")]
    colnames(coords) <- c("x", "y")
    
  } else if (all(c("x", "y") %in% colnames(coords))) {
    
    ## 如果本来就叫 x/y，直接使用
    coords <- coords[, c("x", "y")]
    
  } else if (all(c("col", "row") %in% colnames(coords))) {
    
    ## 有些对象中坐标叫 col/row
    coords <- coords[, c("col", "row")]
    colnames(coords) <- c("x", "y")
    
  } else {
    
    ## 如果都不符合，说明需要手动查看坐标列名
    stop("无法识别空间坐标列名")
  }
  
  ## 对齐 RCTD 结果和空间坐标
  common_spots <- intersect(
    rownames(rctd_mat),
    rownames(coords)
  )
}

## 只保留共同 spot
{rctd_mat <- rctd_mat[common_spots, ]
  coords <- coords[common_spots, ]
  
  ## 检查 RCTD 细胞比例矩阵和空间坐标是否完全对齐，必须返回 TRUE
  stopifnot(
    identical(rownames(rctd_mat), rownames(coords))
  )
}

## 修复细胞类型名字
name_map <- data.frame(
  old_name = colnames(rctd_mat),
  new_name = make.names(
    colnames(rctd_mat),
    unique = TRUE
  )
)

## 替换 rctd_mat 的列名
colnames(rctd_mat) <- name_map$new_name

## 保存新旧名字对应关系，方便后续解释结果时知道哪个变量对应哪个真实细胞类型
write.csv(
  name_map,
  "MISTy_celltype_name_map.csv",
  row.names = FALSE
)

## 构建 MISTy 输入矩阵
{misty_input <- as.data.frame(rctd_mat)
  
  ## 构建 MISTy 多视角 views
  misty_views <- create_initial_view(
    misty_input
  )
}

## 添加近邻视角和远邻视角
{misty_views <- misty_views %>%
    
    add_juxtaview(
      positions = coords[, c("x", "y")],
      neighbor.thr = 50,
      prefix = "juxta"
    ) %>%
    
    
    add_paraview(
      positions = coords[, c("x", "y")],
      l = 100,
      prefix = "para"
    )
}

## 运行 MISTy
misty_result <- run_misty(
  misty_views,
  results.folder = "MISTy_RCTD_result"
)

#保存
save(misty_result,file = "misty_result.rda")
load("misty_result.rda")

## 收集 MISTy 结果
misty_collect <- collect_results(
  "MISTy_RCTD_result"
)

##  绘制不同空间视角贡献图
## 展示每个 target 的模型解释能力，可以判断：哪些细胞类型的空间分布最容易被自身或邻域解释
plot_improvement_stats(
  misty_collect
)

## 展示 predictor 和 target 之间的影响关系
{plot_interaction_heatmap(
  misty_collect,
  view = "intra"
)}

## 自动匹配上皮细胞和巨噬细胞，grep 用于根据关键词搜索列名
ECM_Fib_cell <- grep(
  "PI16_Fib|PI16",
  colnames(rctd_mat),
  value = TRUE
)[1]

Activated_Fib_cell <- grep(
  "Activated_Fib|Activated",
  colnames(rctd_mat),
  value = TRUE
)[1]

## 定义一个通用函数：
## 输入一个细胞类型名，就画该细胞类型在空间中的比例分布
{plot_spatial_celltype <- function(celltype) {
  
  ggplot(
    plot_df,
    aes(
      x = x,
      y = y,
      color = .data[[celltype]]
    )
  ) +
    geom_point(
      size = 1.2,
      alpha = 0.95
    ) +
    scale_color_viridis_c(
      option = "turbo",
      name = celltype
    ) +
    coord_fixed() +
    
    ## 因为空间图像的 y 轴方向通常和 ggplot 默认方向相反
    ## 所以要反转 y 轴，保证切片方向正常
    scale_y_reverse() +
    theme_void(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    ) +
    labs(
      title = celltype
    )
}
}
plot_df <- cbind(
  coords,
  rctd_mat
)

head(plot_df)
## 分别绘制上皮细胞和巨噬细胞
p_PI16_Fib <- plot_spatial_celltype(
  PI16_Fib_cell
)

p_Activated_Fib <- plot_spatial_celltype(
  Activated_Fib_cell
)

## patchwork 拼图
p_final <- p_PI16_Fib + p_Activated_Fib

print(p_final)


##  细胞类型平均比例柱状图
## 计算每种细胞类型在所有 spot 中的平均比例
{cell_prop_df <- data.frame(
  celltype = colnames(rctd_mat),
  mean_prop = colMeans(rctd_mat)
)
  
  ## 绘制平均比例柱状图，可以看哪些细胞类型整体比例更高
  p_bar <- ggplot(
    cell_prop_df,
    aes(
      x = reorder(celltype, mean_prop),
      y = mean_prop
    )
  ) +
    geom_col(
      fill = "grey60"
    ) +
    coord_flip() +
    theme_classic(base_size = 14) +
    labs(
      x = "Cell type",
      y = "Mean proportion",
      title = "Mean RCTD cell-type proportion"
    )
}
p_bar


##  两种细胞类型共定位散点图这里默认取前两个细胞类型作为示例> 
colnames(rctd_mat)
# 载入需要的命名空间：spacexr
{## 实际分析中建议手动指定：
  ## celltype_x <- "Epithelial.cells"
  ## celltype_y <- "Macrophage.cells"
  celltype_x <- colnames(rctd_mat)[14]
  celltype_y <- colnames(rctd_mat)[5]
}

## 绘制两个细胞类型比例之间的散点图
## 如果点呈上升趋势，说明二者可能空间共定位， 如果呈下降趋势，说明可能空间互斥
p_scatter <- ggplot(
  plot_df,
  aes(
    x = .data[[celltype_x]],
    y = .data[[celltype_y]]
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = celltype_x,
    y = celltype_y,
    title = paste0(
      celltype_x,
      " vs ",
      celltype_y
    )
  )
library(ggpubr)
p_scatter+ 
  stat_cor(method = "pearson",   # 或 "spearman"
           label.x.npc = "center", 
           label.y.npc = "top",
           size = 4)
# install.packages("ggpubr")