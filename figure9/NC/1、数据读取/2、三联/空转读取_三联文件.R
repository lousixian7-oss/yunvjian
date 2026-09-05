
library(Seurat)
library(ggplot2)


# 读取基因表达矩阵
counts <- Read10X(data.dir = "counts")
# 创建Seurat对象
stRNA <- CreateSeuratObject(counts = counts, assay = "Spatial", project = "Visium")

# 读取空间图像及坐标信息
{
  image <- Read10X_Image(image.dir = "spatial",filter.matrix = TRUE)
}

# 将图像添加进 Seurat 对象
stRNA[["slice1"]] <- image
DefaultAssay(stRNA) <- "Spatial"

#数据质控可视化
{VlnPlot(stRNA, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1)
  SpatialFeaturePlot(stRNA, features = "nCount_Spatial") + theme(legend.position = "right")
  SpatialFeaturePlot(stRNA, features = "nFeature_Spatial") + theme(legend.position = "right")
}

#归一化与降维分析
{stRNA <- SCTransform(stRNA, assay = "Spatial", verbose = TRUE)
  stRNA <- RunPCA(stRNA, assay = "SCT", verbose = TRUE)
  stRNA <- FindNeighbors(stRNA, dims = 1:15)
  stRNA <- FindClusters(stRNA)
}

#可以修改下面代码，降低分群数目
stRNA <- FindClusters(stRNA, resolution = 0.5)
stRNA <- RunUMAP(stRNA, dims = 1:15)

# 聚类与空间可视化
{DimPlot(stRNA, reduction = "umap", label = TRUE)
  SpatialPlot(stRNA, label = TRUE, label.size = 4)
}

# 基因表达展示
{
  SpatialFeaturePlot(stRNA, features = c("COX4I1"))
  
}
# 保存结果
save(stRNA, file = "yjsl_示例.rda")
