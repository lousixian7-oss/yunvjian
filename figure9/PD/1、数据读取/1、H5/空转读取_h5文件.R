library(Seurat)
library(ggplot2)
library(patchwork)


##读取表达矩阵
## Read10X_h5() 用于读取10X Genomics输出的h5格式表达矩阵
counts <- Read10X_h5(
  filename = "filtered_feature_bc_matrix.h5"
)

##创建Seurat空间对象
{stRNA <- CreateSeuratObject(
  counts = counts,
  assay = "Spatial",
  project = "Visium"
)}

## 读取空间图像和spot坐标
## 当前目录中通常需要包含：
## ├── tissue_positions.csv 或 tissue_positions_list.csv
## ├── scalefactors_json.json
## ├── tissue_lowres_image.png
## └── tissue_hires_image.png
image <- Read10X_Image(
  image.dir = ".",# image.dir = "." 表示空间相关文件在当前工作目录
  filter.matrix = TRUE
)

##将图像信息与Seurat对象中的spot对齐
image <- image[Cells(stRNA)]

## 将空间图像加入Seurat对象
## 这里把image命名为 slice1
## 以后SpatialDimPlot / SpatialFeaturePlot 都会使用这个图像
stRNA[["slice1"]] <- image

## 设置默认assay
## 和单细胞类似，告诉Seurat后续默认使用 Spatial 这个表达矩阵
DefaultAssay(stRNA) <- "Spatial"

##查看原始空间spot分布
## SpatialDimPlot() 将spot显示在组织图像上
## 这里还没有聚类，所以主要用于检查：图像是否成功导入；spot是否正确覆盖在组织区域

{SpatialDimPlot(stRNA)
  }
## 数据质控可视化， VlnPlot() 用小提琴图查看不同spot的质量分布
VlnPlot(
  stRNA,
  features = c("nCount_Spatial", "nFeature_Spatial"),
  pt.size = 0.1
)

## 空间展示每个spot的UMI数量，颜色越深/越高，说明该spot总UMI越多

{SpatialFeaturePlot(
  stRNA,
  features = "nCount_Spatial"
) +
    theme(
      legend.position = "right"
    )
}

## 空间展示每个spot检测到的基因数
{SpatialFeaturePlot(
  stRNA,
  features = "nFeature_Spatial"
) +
    theme(
      legend.position = "right"
    )
  
  
}

## SCTransform标准化
## SCTransform() 是Seurat推荐的标准化方法
## 作用：1消除不同spot测序深度差异，2稳定表达方差，3生成新的 SCT assay
{stRNA <- SCTransform(
  stRNA,
  assay = "Spatial",
  verbose = FALSE
)}

##  PCA降维
## RunPCA() 对SCT标准化后的表达矩阵进行主成分分析
{stRNA <- RunPCA(
  stRNA,
  assay = "SCT",
  verbose = FALSE
)
}

## ElbowPlot() 用来判断选择多少个PC比较合适
ElbowPlot(stRNA)

##  构建邻近图
## FindNeighbors() 根据PCA结果构建spot之间的近邻关系，是聚类前的必要步骤
{stRNA <- FindNeighbors(
  stRNA,
  reduction = "pca",## reduction = "pca"：使用PCA空间
  dims = 1:15## dims = 1:15：使用前15个主成分
)
}

## 聚类
## FindClusters() 根据邻近图对spot进行聚类
stRNA <- FindClusters(
  stRNA,
  resolution = 0.2## resolution 控制聚类颗粒度： resolution 越大，cluster越多，resolution 越小，cluster越少
)

##  UMAP降维
## RunUMAP() 将PCA结果进一步压缩到二维，用于可视化spot在表达空间中的关系
{
  stRNA <- RunUMAP(
    stRNA,
    reduction = "pca",
    dims = 1:15
  )
}

##  UMAP聚类可视化
## DimPlot() 在UMAP图上显示不同cluster
DimPlot(
  stRNA,
  reduction = "umap",
  label = TRUE## label = TRUE 会在图中标注cluster编号
)

## 空间聚类可视化
## SpatialDimPlot() 在组织切片上显示cluster分布
## 这个图非常重要：它能显示不同转录组cluster是否具有空间结构
{SpatialDimPlot(
  stRNA,
  label = TRUE,
  label.size = 4
)
}
## 单基因空间表达展示
## SpatialFeaturePlot() 可显示某个基因在组织中的空间表达
## 这里展示 MMP14， MMP14 通常与基质重塑、侵袭、迁移有关
SpatialFeaturePlot(stRNA,
                   features = c("MMP14"))

## 保存分析后的stRNA对象
save(stRNA,file = "h5_示例数据.rda")