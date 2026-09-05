library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
## 读取前面保存好的空间转录组Seurat对象
load("yjsl_示例.rda")
##  查看原始聚类结果


#可以修改下面代码，降低分群数目
stRNA <- FindClusters(stRNA, resolution = 0.6)
stRNA <- RunUMAP(stRNA, dims = 1:15)

# 聚类与空间可视化
{DimPlot(stRNA, reduction = "umap", label = TRUE)
  SpatialPlot(stRNA, label = TRUE, label.size = 4)
}

{p1 <- SpatialDimPlot(
  stRNA,
  label = TRUE,
  label.size = 5
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
all_markers<-read.csv("all_markers.csv")
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

# 1. 上层（3,4,9）→ 完全分化的基底上层上皮（包含角化/非角化）
stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("3", "4", "9")] <- "Suprabasal epithelial"

# 2. 中层炎症（8,11） + 下层炎症/重塑（6）→ 炎症上皮
stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c( "8", "11")] <- "Inflamed epithelium"

# 3. 中层（1,5,0）→ 基底/增殖上皮（位于中间层）
stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("0", "1", "5")] <- "Basal epithelium"

# 4. 下层（2）→ 网状基质（成纤维细胞）
stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("2","10")] <- "Reticular stroma"

# 5. 下层（7）→ 免疫-基质（B/浆细胞浸润）
stRNA@meta.data$Region[stRNA@meta.data$seurat_clusters %in% c("6","7")] <- "Immune stroma"


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
  pt.size.factor = 2,
  cols = c(
    "Basal epithelium"           = "#EF9B94",  # 暖鲑鱼粉 (原 #FDDED7 微调)
    "Suprabasal epithelial"      = "#EAA369",  # 暖杏橙 (原 #F5BE8F 微调)
    "Inflamed epithelium"        = "#D35D5D",  # 醒目砖红 (汲取 #e3632d 与 #CC79A7 调优)
    "Reticular stroma"          = "#72B3AD",  # 冷蓝绿 (原 #C1E0DB 调深)
    "Immune stroma"             = "#6B8DAE"  # 柔和淡紫 (原 #A28CC2 微调)
    
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
    "Basal epithelium"           = "#EF9B94",  # 暖鲑鱼粉 (原 #FDDED7 微调)
    "Suprabasal epithelial"      = "#EAA369",  # 暖杏橙 (原 #F5BE8F 微调)
    "Inflamed epithelium"        = "#D35D5D",  # 醒目砖红 (汲取 #e3632d 与 #CC79A7 调优)
    "Reticular stroma"          = "#72B3AD",  # 冷蓝绿 (原 #C1E0DB 调深)
    "Immune stroma"             = "#6B8DAE"  # 柔和淡紫 (原 #A28CC2 微调)
    
  )
)
}

## 单个基因的空间表达展示；展示某个基因在组织切片上的空间表达分布
## 这里以MMP14为例
{SpatialFeaturePlot(
  stRNA,
  features = c("COX4I1")
)
}
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
## 每个Region筛选Top10 marker基因
{top10_region_markers <- region_markers %>%
    group_by(cluster) %>%## 按Region分组
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%## 在每个Region内部按照avg_log2FC从高到低排序
    slice_head(n = 5) %>%## 每个Region取前10个marker基因
    ungroup()
}

## unique()用于去重，避免不同Region中出现重复基因
top10_genes <- unique(top10_region_markers$gene)

##  绘制Region marker热图
DefaultAssay(stRNA) <- "SCT"

{top10_genes <- intersect(top10_genes, rownames(stRNA))
  
  stRNA <- ScaleData(
    stRNA,
    assay = "SCT",
    features = top10_genes
  )
  
  stRNA$Region <- factor(
    stRNA$Region,
    levels = c( "Basal epithelium" ,  
                "Suprabasal epithelial" ,  
                "Inflamed epithelium", 
                "Reticular stroma",  
                "Immune stroma"  
    )
  )
  
  Idents(stRNA) <- "Region"
}
##绘图
{p_heatmap <- DoHeatmap(
  stRNA,
  features = top10_genes,
  group.by = "Region",
  assay = "SCT",
  size = 3.5,
  angle = 0,
  draw.lines = FALSE,
  group.bar = TRUE,
  group.colors = c(
    "Basal epithelium"           = "#EF9B94",  # 暖鲑鱼粉 (原 #FDDED7 微调)
    "Suprabasal epithelial"      = "#EAA369",  # 暖杏橙 (原 #F5BE8F 微调)
    "Immune epithelium"        = "#D35D5D",  # 醒目砖红 (汲取 #e3632d 与 #CC79A7 调优)
    "Reticular stromal"          = "#72B3AD",  # 冷蓝绿 (原 #C1E0DB 调深)
    "Immune stromal"             = "#6B8DAE"  # 柔和淡紫 (原 #A28CC2 微调)
   
  )
) +
    scale_fill_gradientn(
      colors = c("#313695", "#F7F7F7", "#A50026"),
      limits = c(-2, 2),
      oob = scales::squish,
      name = "Scaled\nexpression"
    ) +
    ggtitle("Spatial region-specific marker genes") +
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
  file = "牙龈NC_注释.rda"
)
load("牙龈NC_注释.rda")

NC <- stRNA

load("牙周炎_注释.rda")

PD <- stRNA
genes_show <- c(
  "C3",
  "AGT",
  "PI16",
  "COL1A1",
  "FN1",
  "POSTN",
  "CD44"
)
SpatialFeaturePlot(
  PD,
  features=genes_show,
  pt.size.factor=5,
  ncol=7
)
SpatialFeaturePlot(
  NC,
  features=genes_show,
  pt.size.factor=5,
  ncol=7
)

DefaultAssay(PD) <- "Spatial"
DefaultAssay(NC) <- "Spatial"


PD <- NormalizeData(
  PD,
  assay="Spatial"
)

NC <- NormalizeData(
  NC,
  assay="Spatial"
)
library(dplyr)
library(tidyr)


get_spatial_bubble <- function(obj, group_name, genes){
  
  expr <- GetAssayData(
    obj,
    assay="Spatial",
    layer="data"
  )
  
  expr <- expr[genes, , drop=FALSE]
  
  
  df <- data.frame(
    gene = rownames(expr),
    avg_exp = rowMeans(expr),
    pct_exp = rowMeans(expr > 0) * 100
  )
  
  
  df$Group <- group_name
  
  return(df)
}


PD_df <- get_spatial_bubble(
  PD,
  "PD",
  genes_show
)


NC_df <- get_spatial_bubble(
  NC,
  "NC",
  genes_show
)


bubble_df <- rbind(
  PD_df,
  NC_df
)
library(ggplot2)


p_bubble <- ggplot(
  bubble_df,
  aes(
    x=gene,
    y=Group
  )
)+
  geom_point(
    aes(
      size=pct_exp,
      color=avg_exp
    )
  )+
  scale_color_viridis_c(
    option="plasma"
  )+
  scale_size(
    range=c(2,10)
  )+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle=45,
      hjust=1,
      size=12
    ),
    axis.text.y = element_text(
      size=12
    ),
    axis.title=NULL
  )+
  labs(
    color="Average expression",
    size="% spots expressing"
  )


p_bubble