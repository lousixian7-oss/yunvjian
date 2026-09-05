library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
## 读取前面保存好的空间转录组Seurat对象
load("牙周炎.rda")
##  查看原始聚类结果

#可以修改下面代码，降低分群数目
stRNA <- FindClusters(stRNA, resolution = 0.8)
stRNA <- RunUMAP(stRNA, dims = 1:15)

# 聚类与空间可视化
{DimPlot(stRNA, reduction = "umap", label = TRUE)
  SpatialPlot(stRNA, label = TRUE, label.size = 4)
}

{p1 <- SpatialDimPlot(
  stRNA,
  label = TRUE,
  label.size = 3
)
}

p1
## DimPlot：在UMAP降维图上展示Seurat聚类结果
{p2 <- DimPlot(
  stRNA,
  reduction = "umap",
  label = TRUE
)
}
p2
## 左边通常是空间图，右边是UMAP图
p1 + p2
# 确定当前使用的聚类结果（假设存在一个名为 "seurat_clusters" 的metadata）
Idents(stRNA) <- "seurat_clusters"

# 找出所有簇的标记基因（只保留阳性标记，logFC阈值0.25，最小表达比例0.1）
all_markers <- FindAllMarkers(stRNA, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(all_markers, file = "all_markers.csv", row.names = FALSE)
# 查看每个簇的前5个标记基因
top_markers <- all_markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)

# 打印出来
print(top_markers,n=215)


## 注释方法可以很灵活：可以根据空间位置、marker基因、病理区域、聚类结果等综合判断
## 这里是把已有的seurat_clusters重新归并成三个空间区域：Type1 / Type2 / Type3
## 先新建一个Region列，初始值设为NA
stRNA@meta.data$Region <- NA

# 初始化 Region 列（如果不存在）
if (!"Region" %in% colnames(stRNA@meta.data)) {
  stRNA@meta.data$Region <- NA
}


stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("2")] <- "Suprabasal keratinized"

stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("6")] <- "Suprabasal granular"

stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("3")] <- "Suprabasal spinous"

stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("5")] <- "Basal epithelium"

stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("4")] <- "Junctional epithelium"

stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("1")] <- "Reticular stroma"

stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("0")] <- "Immune stroma"


# （可选）检查是否所有簇都被覆盖
unique(stRNA@meta.data$Region[!is.na(stRNA@meta.data$Region)])
## 检查每个区域有多少spot
table(stRNA$Region)


## 空间区域可视化
## SpatialDimPlot：在组织切片上展示Region区域分布

{p_region <- SpatialDimPlot(
  stRNA,
  group.by = "Region",
  label = F,
  label.size = 3,
  repel = TRUE,
  pt.size.factor = 4.5,
  cols = c(
    "Suprabasal keratinized"  = "#F7D9B0",  # 浅米杏（表层角化）
    "Suprabasal granular"     = "#EAA369",  # 暖杏橙（颗粒层）
    "Suprabasal spinous"      = "#D99B5A",  # 深杏/橘棕（棘层）
    "Basal epithelium"        = "#EF9B94",  # 鲑鱼粉（基底层，与健康一致）
    "Junctional epithelium"   = "#D35D5D",  # 砖红（结合上皮，炎症特征）
    "Immune stroma"           = "#6B8DAE",  # 冷灰蓝（免疫基质，与健康一致）
    "Reticular stroma"        = "#72B3AD"   # 冷蓝绿（网状基质，与健康一致）
    
  )
) +
    ggtitle("") +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "right"
    )
}
## 显示空间区域注释图
p_region

## UMAP上展示Region注释结果
{DimPlot(
  stRNA,
  reduction = "umap",
  group.by = "Region",
  label = F,
  cols = c(
    "Suprabasal keratinized"  = "#F7D9B0",  # 浅米杏（表层角化）
    "Suprabasal granular"     = "#EAA369",  # 暖杏橙（颗粒层）
    "Suprabasal spinous"      = "#D99B5A",  # 深杏/橘棕（棘层）
    "Basal epithelium"        = "#EF9B94",  # 鲑鱼粉（基底层，与健康一致）
    "Junctional epithelium"   = "#D35D5D",  # 砖红（结合上皮，炎症特征）
    "Immune stroma"           = "#6B8DAE",  # 冷灰蓝（免疫基质，与健康一致）
    "Reticular stroma"        = "#72B3AD"   # 冷蓝绿（网状基质，与健康一致）
  )
)
}

## 单个基因的空间表达展示；展示某个基因在组织切片上的空间表达分布

DefaultAssay(stRNA) <- "Spatial"

stRNA <- NormalizeData(
  stRNA,
  assay = "Spatial",
  normalization.method = "LogNormalize"
)

Layers(stRNA[["Spatial"]])
SpatialFeaturePlot(
  stRNA,
  features = c(
    "PI16",
    "AGT",
    "C3",
    "COL1A1",
    "FN1",
    "POSTN"
  ),
  pt.size.factor = 5,
  max.cutoff = "q95",
  ncol = 6          # 新增，六张图排成一行
)
SpatialFeaturePlot(
  stRNA,
  c(
    "BRCA1",
    "CYP3A4",
    "CXCR4",
    "FOS"
  ),
  pt.size.factor = 5,
  max.cutoff = 'q95'
) 


#+
#   scale_fill_gradientn(colors = c("lightgrey", "yellow"), limits = c(0, 6))

## 差异基因分析准备
library(dplyr)
## 设置Seurat身份标签，后续FindAllMarkers会按照Region分组寻找marker基因
Idents(stRNA) <- "Region"

## 使用SCT assay做差异分析
DefaultAssay(stRNA) <- "SCT"

## 寻找Type1 / Type2 / Type3各自的marker基因
## FindAllMarkers：对每个Region分别寻找相对于其他所有Region高表达的基因

{region_markers <- FindAllMarkers(
  stRNA,
  assay = "SCT",
  only.pos = TRUE,## 只保留在该区域上调的marker基因
  min.pct = 0.25,## 该基因至少在某一组25%的spot中表达
  logfc.threshold = 0.25## 平均logFC至少大于0.25才纳入结果
)
}
## 保存Region marker基因
write.csv( region_markers,file = "Region_all_markers.csv",row.names = FALSE)
region_markers<-read.csv("Region_all_markers.csv")
## 每个Region筛选Top5 marker基因
{top5_region_markers <- region_markers %>%
    group_by(cluster) %>%## 按Region分组
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%## 在每个Region内部按照avg_log2FC从高到低排序
    slice_head(n = 5) %>%## 每个Region取前5个marker基因
    ungroup()
}

## unique()用于去重，避免不同Region中出现重复基因
top5_genes <- unique(top5_region_markers$gene)

##  绘制Region marker热图
DefaultAssay(stRNA) <- "SCT"

{top5_genes <- intersect(top5_genes, rownames(stRNA))
  
  stRNA <- ScaleData(
    stRNA,
    assay = "SCT",
    features = top5_genes
  )
  
  stRNA$Region <- factor(
    stRNA$Region,
    levels = c("Suprabasal keratinized",   # 表层角化
               "Suprabasal granular",      # 颗粒层
               "Suprabasal spinous",       # 棘层
               "Basal epithelium",         # 基底层
               "Junctional epithelium",    # 结合上皮（牙周炎特有）
               "Immune stroma",            # 免疫基质（浅层基质）
               "Reticular stroma"          # 网状基质（深层基质）
               )
  )
  
  Idents(stRNA) <- "Region"
}
##绘图
{p_heatmap <- DoHeatmap(
  stRNA,
  features = top5_genes,
  group.by = "Region",
  assay = "SCT",
  size = 3.5,
  angle = 45,
  draw.lines = FALSE,
  group.bar = TRUE,
  group.colors = c(
    "Suprabasal keratinized"  = "#F7D9B0",  # 浅米杏（表层角化）
    "Suprabasal granular"     = "#EAA369",  # 暖杏橙（颗粒层）
    "Suprabasal spinous"      = "#D99B5A",  # 深杏/橘棕（棘层）
    "Basal epithelium"        = "#EF9B94",  # 鲑鱼粉（基底层，与健康一致）
    "Junctional epithelium"   = "#D35D5D",  # 砖红（结合上皮，炎症特征）
    "Immune stroma"           = "#6B8DAE",  # 冷灰蓝（免疫基质，与健康一致）
    "Reticular stroma"        = "#72B3AD"   # 冷蓝绿（网状基质，与健康一致）
  )
) +
    scale_fill_gradientn(
      colors = c("#313695", "#F7F7F7", "#A50026"),
      limits = c(-2, 2),
      oob = scales::squish,
      name = "Scaled\nexpression"
    ) +
    ggtitle("") +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9, face = "italic", color = "black"),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "right"
    )
}
p_heatmap

## 保存对象后，后续分析可以直接加载：
save(stRNA,
  file = "牙周炎_注释.rda"
)

load("牙周炎_注释.rda")