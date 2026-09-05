library(Seurat)
library(spacexr)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(ggpubr)
library(viridis)
library(scales)


## 读取数据
## 单细胞参考数据，要求：
## 1）是Seurat对象
## 2）已经完成细胞类型注释
## 3）meta.data中有 cellType 这一列

yjsl<-readRDS("yjsl_Type.rds")

# 4. 检查结果
table(yjsl$cellType)

load("牙周炎_注释.rda")
##  提取单细胞表达矩阵
## 设置默认assay为RNA
DefaultAssay(yjsl) <- "RNA"

## 提取单细胞原始counts矩阵，行是基因，列是单细胞barcode
{sc_counts <- GetAssayData(
  yjsl,
  assay = "RNA",
  layer = "counts"
)
}

##  提取单细胞细胞类型标签，cellType 是每个单细胞的注释结果
## 例如：
## Epithelial cells
## Macrophage cells
## T cells
## Fibroblasts
cell_types <- yjsl@meta.data$cellType

## 给cell_types加上细胞barcode名称
{names(cell_types) <- colnames(yjsl)
}
## 查看每种细胞类型数量
## 如果某种细胞太少，RCTD参考会不稳定
table(cell_types)


# 构建 Reference 对象
reference <- Reference(
  counts = sc_counts,
  cell_types = factor(cell_types)
)


## 提取空间转录组表达矩阵
## 设置默认assay为Spatial
# DefaultAssay(stRNA) <- "Spatial"
# 
# ## 提取空间转录组counts，行是基因，列是空间spot barcode
# st_counts <- GetAssayData(
#   stRNA,assay = "Spatial",layer = "counts"
# )
# 
# ## 提取空间坐标， GetTissueCoordinates() 提取每个spot在组织切片上的空间坐标
# {coords <- GetTissueCoordinates(stRNA)
# }
# 
# {## 不同Seurat版本坐标列名可能不同
#   ## 如果没有x/y，需要根据实际情况改成imagecol/imagerow等
#   coords <- coords[, c("x", "y")]
#   
#   ## 让坐标顺序与表达矩阵spot顺序完全一致
#   ## RCTD要求：st_counts的列顺序 = coords的行顺序
#   coords <- coords[colnames(st_counts), ]
# }
# 
# ## 检查表达矩阵和坐标是否匹配，必须返回 TRUE
# ## 如果是 FALSE，说明spot顺序或名称不一致，后面不能继续
# all(colnames(st_counts) == rownames(coords))
# 
# 
# ## 构建SpatialRNA对象
# {puck <- SpatialRNA(
#   coords = coords,
#   counts = st_counts
# )
# }
DefaultAssay(stRNA) <- "Spatial"
st_counts <- GetAssayData(stRNA, assay = "Spatial", layer = "counts")

# 4.1 从图像槽中提取正确的坐标（修复行名为NA的问题）
# 尝试从 boundaries 中获取 centroids
coord_df <- stRNA@images$slice1@boundaries$centroids@coords
barcodes <- stRNA@images$slice1@boundaries$centroids@cells
rownames(coord_df) <- barcodes
coords <- coord_df[, c("x", "y")]   # 列名可能不同，实际通常为 "x", "y"

# 4.2 只保留组织内的 spots（即坐标表中存在的spot）
common <- intersect(colnames(st_counts), rownames(coords))
cat("保留的组织内spot数量:", length(common), "\n")

st_counts <- st_counts[, common, drop = FALSE]
coords <- coords[common, , drop = FALSE]

# 4.3 确保顺序完全一致
stopifnot(all(colnames(st_counts) == rownames(coords)))

# 4.4 将坐标转换为 data.frame（SpatialRNA 所需格式）
coords <- as.data.frame(coords)
colnames(coords) <- c("x", "y")

# 4.5 创建 SpatialRNA 对象
nUMI <- colSums(st_counts)
puck <- SpatialRNA(coords = coords, counts = st_counts, nUMI = nUMI)
## 创建RCTD对象
# create.RCTD() 将空间对象和单细胞参考对象组合起来

myRCTD <- create.RCTD(
  puck,
  reference,
  max_cores = 12## 使用CPU核心数，加快计算速度
)

## 运行RCTD反卷积（时间偏长）
## 常用模式：
## "doublet"：适合10X Visium，推荐
## "full"：允许更多细胞类型混合，但更慢
## "singlet"：假设每个spot只有一种细胞类型
{myRCTD <- run.RCTD(myRCTD,
                    doublet_mode = "doublet"
)
}

## 保存RCTD对象
save(myRCTD,file = "myRCTD合并.rda")
load("myRCTD合并.rda")
## 提取RCTD主要预测结果
results <- myRCTD@results$results_df

head(results)

##将RCTD主导细胞类型加入Seurat空间对象
{stRNA$RCTD_celltype <- results$first_type
  
  ##提取RCTD细胞比例矩阵
  weights <- myRCTD@results$weights
}
##  将特定细胞类型比例加入Seurat对象
## 注意：这里的列名必须和 weights 中的细胞类型名称完全一致
## 先运行 colnames(weights) 查看真实名称
colnames(weights)

##换成自己感兴趣的细胞亚群名
## 添加上皮细胞比例
stRNA$Endothelial <- weights[, "Epithelial cells"]
stRNA$Tcell <- weights[, "T cells"]
## 添加巨噬细胞比例
stRNA$Macrophage <- weights[, "Macrophages"]

## 绘制上皮细胞空间比例图

## SpatialFeaturePlot：
## 在组织空间上展示某个数值变量
## 这里展示每个spot的上皮细胞比例
{p_epi <- SpatialFeaturePlot(
  stRNA,
  features = "Endothelial",
  alpha = c(0.01, 1),pt.size.factor = 3
) +
    scale_fill_gradientn(
      colors = c(
        "#F7FBFF",
        "#FFFF00"
        
      ) ,trans = "sqrt"
    ) 
    theme_classic()
}
p_epi
p_fib<- SpatialFeaturePlot(
  stRNA,
  features = "Tcell",
  alpha = c(0.01, 1),pt.size.factor = 3
) +
  scale_fill_gradientn(
    colors = c("#F7FBFF",
               "#FFFF00")
  ) +
  theme_classic()

print(p_fib)
## 绘制巨噬细胞空间比例图
{p_mac <- SpatialFeaturePlot(
  stRNA,
  features = "Macrophage",
  alpha = c(0.01, 1),pt.size.factor = 3
) +
    scale_fill_gradientn(
      colors = c(
        "#F7FBFF",
        "#FFFF00"  # 深红，比例高
      )
    ) +
    theme_classic()
}

p_mac

## 不同Region中细胞比例比较
{library(ggpubr)
  
  ## 提取meta.data用于统计绘图
  df <- stRNA@meta.data
}
## 比较不同Region中的上皮细胞比例

ggboxplot(
  df,
  x = "Region",
  y = "Endothelial",
  fill = "Region"
) +
  stat_compare_means() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
  )



##更高级的表现手法
weights <- as.data.frame(myRCTD@results$weights)

# 归一化为比例
weights <- weights / rowSums(weights)
weights[is.na(weights)] <- 0

# 对齐spot
common_spots <- intersect(rownames(weights), rownames(stRNA@meta.data))
weights <- weights[common_spots, , drop = FALSE]

meta_df <- stRNA@meta.data[common_spots, , drop = FALSE]
meta_df$spot_id <- rownames(meta_df)

## 不同Region平均细胞组成堆叠柱状图

{plot_df <- cbind(
  meta_df[, c("spot_id", "Region")],
  weights
)
  
  plot_long <- plot_df %>%
    pivot_longer(
      cols = -c(spot_id, Region),
      names_to = "Celltype",
      values_to = "Proportion"
    )
  
  avg_celltype <- plot_long %>%
    group_by(Region, Celltype) %>%
    summarise(
      Mean_proportion = mean(Proportion, na.rm = TRUE),
      .groups = "drop"
    )
  
  p_stack <- ggplot(
    avg_celltype,
    aes(x = Region, y = Mean_proportion, fill = Celltype)
  ) +
    geom_col(width = 0.75, color = "white", linewidth = 0.2) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 17)
    ) +
    labs(
      title = "Cell-type composition across spatial regions",
      x = NULL,
      y = "Mean proportion",
      fill = "Cell type"
    )
}

p_stack

## Region × Celltype 平均比例热图
avg_mat <- avg_celltype %>%
  pivot_wider(
    names_from = Region,
    values_from = Mean_proportion
  ) %>%
  as.data.frame()

rownames(avg_mat) <- avg_mat$Celltype
avg_mat$Celltype <- NULL
avg_mat <- as.matrix(avg_mat)

pheatmap(
  avg_mat,
  scale = "row",
  color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
  border_color = NA,
  clustering_method = "ward.D2",
  fontsize_row = 10,
  fontsize_col = 12,
  main = "cell-type enrichment"
)

##所有细胞类型在Region中的小提琴图

p_violin_all <- ggplot(
  plot_long,
  aes(x = Region, y = Proportion, fill = Region)
) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ Celltype, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = c(
      "Suprabasal keratinized"  = "#F7D9B0",  # 浅米杏（表层角化）
      "Suprabasal granular"     = "#EAA369",  # 暖杏橙（颗粒层）
      "Suprabasal spinous"      = "#D99B5A",  # 深杏/橘棕（棘层）
      "Basal epithelium"        = "#EF9B94",  # 鲑鱼粉（基底层，与健康一致）
      "Junctional epithelium"   = "#D35D5D",  # 砖红（结合上皮，炎症特征）
      "Immune stroma"           = "#6B8DAE",  # 冷灰蓝（免疫基质，与健康一致）
      "Reticular stroma"        = "#72B3AD"   # 冷蓝绿（网状基质，与健康一致）
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  labs(
    title = "Cell-type proportions across spatial regions",
    x = NULL,
    y = "RCTD proportion"
  )

p_violin_all

