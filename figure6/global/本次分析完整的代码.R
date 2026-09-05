
library(Seurat)#单细胞数据分析
library(dplyr)
library(patchwork)#将多个ggplot图像组合成一个图
library(ggplot2)
library(SingleR)#对单细胞数据进行自动注释
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(GSEABase)
library(harmony)
library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)

library(UCell)
library(AUCell)
library(decoupleR)
library(singscore)

library(irGSEA)##安装比较麻烦
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)

library(rhdf5)
################读取并构建seurat对象################
################ 读取并构建 Seurat 对象（BD Rhapsody 格式）################
{
  
  data_dir <- paste0(getwd(), "/data")       # 数据根目录
  samples <- list.files(data_dir)            # 所有样本文件夹名
  seurat_list <- list()                      # 用于存放每个样本的 Seurat 对象
  
  for (s in samples) {
    sample_path <- file.path(data_dir, s)
    
    # 找到该样本文件夹下的 counts 文件和 cellname 文件
    count_file <- list.files(sample_path, pattern = "counts.tsv.gz$", full.names = TRUE)
    cell_file  <- list.files(sample_path, pattern = "cellname.list.txt.gz$", full.names = TRUE)
    
    if (length(count_file) == 0) {
      warning(paste("未找到 counts 文件，跳过样本：", s))
      next
    }
    
    # 1. 读取表达矩阵（行为基因，列为细胞）
    counts <- read.table(count_file[1], 
                         header = TRUE, 
                         row.names = 1, 
                         check.names = FALSE)
    
    # 2. 如果有 cellname 文件，用其重命名细胞（确保列名准确）
    if (length(cell_file) > 0) {
      cell_names <- read.table(cell_file[1], header = FALSE, stringsAsFactors = FALSE)$V1
      if (ncol(counts) == length(cell_names)) {
        colnames(counts) <- cell_names
      } else {
        warning(paste("细胞数量与 cellname 文件不匹配，样本：", s))
      }
    }
    
    # 3. 创建 Seurat 对象（先不过滤，后面统一质控）
    seurat_obj <- CreateSeuratObject(counts = counts,
                                     project = s,
                                     min.cells = 0,      # 暂不过滤基因
                                     min.features = 0)   # 暂不过滤细胞
    
    # 4. 添加样本信息（便于后续区分）
    seurat_obj$Type <- s
    seurat_obj$sample <- s
    
    seurat_list[[s]] <- seurat_obj
  }
  
  # 5. 合并所有样本（如果多于 1 个）
  if (length(seurat_list) > 1) {
    yjsl <- merge(seurat_list[[1]], 
                  y = seurat_list[-1], 
                  add.cell.ids = names(seurat_list))
  } else {
    yjsl <- seurat_list[[1]]
  }
  
  # 6. 查看合并后的结果
  print(yjsl)
}


# 定义需要保留的样本（9个，不含治疗后）
keep_samples <- c("GSM5220920", "GSM5220921", "GSM5220922", "GSM5220923",
                  "GSM5220924", "GSM5220925", "GSM5220926", "GSM5220927",
                  "GSM5220928")

# 用 Type 列过滤
yjsl <- subset(yjsl, subset = Type %in% keep_samples)

# 验证结果
ncol(yjsl)  # 应该是 45191
table(yjsl$Type)
# 根据Type组的样本名字定义新的group组
yjsl$group <- ifelse(yjsl$Type %in% c("GSM5220920", "GSM5220921","GSM5220922","GSM5220923"), "NC", 
                     ifelse(yjsl$Type %in% c("GSM5220925", "GSM5220926","GSM5220927","GSM5220928","GSM5220924"), "PD", NA))

# 查看前几行数据，确认分组信息
head(yjsl@meta.data)
table(yjsl$group)

data2_dir <-  paste0(getwd(), "/data2")  # 路径根据实际调整
samples2 <- list.files(data2_dir)

seurat_list2 <- list()

for (s in samples2) {
  sample_path <- file.path(data2_dir, s)
  # 读取10x格式
  counts <- Read10X(data.dir = sample_path)
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = counts,
                                   project = s,
                                   min.cells = 0,
                                   min.features = 0)
  # 添加样本信息
  seurat_obj$Type <- s
  seurat_obj$sample <- s
  seurat_list2[[s]] <- seurat_obj
}

# 合并所有data2样本
if (length(seurat_list2) > 1) {
  yjsl2 <- merge(seurat_list2[[1]], y = seurat_list2[-1], add.cell.ids = names(seurat_list2))
} else {
  yjsl2 <- seurat_list2[[1]]
}


library(rhdf5)

data3_dir <-  paste0(getwd(), "/data3") 
h5_files <- list.files(data3_dir, pattern = "\\.h5$", full.names = TRUE)

seurat_list3 <- list()

for (file in h5_files) {
  # 提取样本名（如 GSM4600896）
  sample_name <- gsub("_filtered_feature_bc.*$", "", basename(file))
  
  # 读取 h5 文件（注意使用 Read10X_h5）
  counts <- Read10X_h5(filename = file)
  seurat_obj <- CreateSeuratObject(counts = counts,
                                   project = sample_name,
                                   min.cells = 0,
                                   min.features = 0)
  seurat_obj$Type <- sample_name
  seurat_obj$sample <- sample_name
  seurat_list3[[sample_name]] <- seurat_obj
}

# 合并data3样本
if (length(seurat_list3) > 1) {
  yjsl3 <- merge(seurat_list3[[1]], y = seurat_list3[-1], add.cell.ids = names(seurat_list3))
} else {
  yjsl3 <- seurat_list3[[1]]
}


# 首先合并 data2 和 data3（如果两者都有）
if (exists("yjsl2") && exists("yjsl3")) {
  yjsl_all <- merge(yjsl2, y = yjsl3, add.cell.ids = c("data2", "data3"))
} else if (exists("yjsl2")) {
  yjsl_all <- yjsl2
} else if (exists("yjsl3")) {
  yjsl_all <- yjsl3
} else {
  stop("没有可合并的数据")
}

# 再将原来的 yjsl（data1）合并进来
if (exists("yjsl")) {
  yjsl_final <- merge(yjsl, y = yjsl_all, add.cell.ids = c("data1", "others"))
} else {
  yjsl_final <- yjsl_all
}

# 查看合并后的细胞数
print(yjsl_final)


library(dplyr)

yjsl_final$group <- case_when(
  # data1
  yjsl_final$Type %in% c("GSM5220920", "GSM5220921", "GSM5220922", "GSM5220923") ~ "NC",
  yjsl_final$Type %in% c("GSM5220925", "GSM5220926", "GSM5220927", "GSM5220928", "GSM5220924") ~ "PD",
  # data2
  yjsl_final$Type %in% c("GSM5005048", "GSM5005049", "GSM5005050", "GSM5005051", "GSM5005052","GSM5005053", "GSM5005054", "GSM5005055", "GSM5005056", "GSM5005057","GSM5177039","GSM5177040","GSM5177041") ~ "NC",
  yjsl_final$Type %in% c("GSM5005058", "GSM5005059", "GSM5005060", "GSM5005061","GSM5005062", "GSM5177042", "GSM5177043", "GSM5177044") ~ "PD",
  # data3
  yjsl_final$Type %in% c("GSM4600896", "GSM4600897") ~ "NC",
  yjsl_final$Type %in% c("GSM4600898", "GSM4600899") ~ "PD",
  TRUE ~ NA_character_
)

# 检查分组分布
table(yjsl_final$group)


# 查看过滤后的细胞数
ncol(yjsl_final)
yjsl<-yjsl_final


{
################质控可视化################
#计算线粒体基因的百分比
yjsl[["percent.mt"]] <- PercentageFeatureSet(yjsl, pattern = "^MT-")
#计算核糖体基因的百分比
yjsl[["percent.rb"]] <- PercentageFeatureSet(yjsl, pattern = "^RP")
#绘制小提琴图，展示每个细胞中的基因数（nFeature_RNA）、检测到的分子数（nCount_RNA）以及线粒体基因的百分比（percent.mt）
VlnPlot(yjsl, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
        ncol = 3,pt.size = 0)#

##nCount表示每个细胞中检测到的分子数
##nFeature表示每个细胞中检测到的基因数
###上述两者，若有异常的离群细胞则说明该细胞可能为双细胞或多细胞（即，质控目的为剔除异常值）
## 绘制散点图，比较不同特征之间的关系
plot2 <- FeatureScatter(yjsl, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot2
plot3 <- FeatureScatter(yjsl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
plot3

plot2 + plot3
###过滤基因
sub1 <- yjsl$nCount_RNA >= 1000 ## 过滤掉低于1000的细胞
sub2 <- yjsl$nFeature_RNA >= 200 & yjsl$nFeature_RNA <= 8000# 过滤掉基因数小于200或大于8000的细胞
sub3 <- yjsl$percent.mt <= 30# 过滤掉线粒体基因百分比大于20%的细胞
sub4<-yjsl$percent.rb<= 20 # 过滤掉核糖体基因百分比大于20%的细胞
sub <-sub1& sub2 & sub3 # 合并过滤条件
yjsl <- yjsl[, sub]# 应用过滤条件，保留符合条件的细胞
}
#重新绘制
VlnPlot(yjsl, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
        ncol = 3)

table(yjsl$Type)

#######使用LogNormalize方法#######
yjsl <- NormalizeData(yjsl, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

#######鉴定高变基因#######

yjsl <- FindVariableFeatures(yjsl, selection.method = "vst", nfeatures = 2000)

#提取前15的高变基因
top15 <- head(VariableFeatures(yjsl), 15)
top15

##展示高变基因
plot1 <- VariableFeaturePlot(yjsl, pt.size = 2, raster = T)
plot1
plot2 <- LabelPoints(plot = plot1, points = top15, 
                     xnudge = 0,
                     ynudge = 0)
plot2 


######数据Scaling，降维前的必要准备
{
  all.genes <- rownames(yjsl)
  yjsl <- ScaleData(yjsl, features = VariableFeatures(yjsl))
}
####正式降维###
yjsl<- RunPCA(yjsl, features = VariableFeatures(object =yjsl)) 
#数据可视化方法
#VizDimReduction
VizDimLoadings(yjsl, dims = 1:2, reduction = "pca")

#DimHeatmap
DimHeatmap(yjsl, dims = 1:15, cells = 500, balanced = TRUE)
#ElbowPlot
ElbowPlot(yjsl, ndims=20, reduction="pca") 

########harmony去除批次效应########################
yjsl <- RunHarmony(yjsl, group.by.vars = "dataset")

########确定维度、细胞分群########
#选择PC
PC=1:15
###根据上一步中的合适维度对细胞进行分群
yjsl <- FindNeighbors(yjsl, dims = PC, reduction = "harmony")
#resolution可以设0.1-1之间，值越高，亚群数目越多，常规0.5
yjsl <- FindClusters(yjsl, resolution =seq(0.2,0.8,0.1))
clustree(yjsl)


##!!!记住这个维度

yjsl <- FindClusters(yjsl, resolution =0.3)
###使用UMAP及进行非线性降维
yjsl <- RunUMAP(yjsl, dims = PC, reduction = "harmony")
#数据可视化
DimPlot(yjsl, reduction = "umap", label = T,repel = T)

DimPlot(yjsl, reduction = "umap", label = T,repel = T,split.by = 'group')

# ###使用tSNE进行非线性降维
# yjsl = RunTSNE(yjsl, dims = PC)
# embed_tsne <- Embeddings(yjsl, 'tsne')
# DimPlot(yjsl, reduction = "tsne") 
# DimPlot(yjsl, reduction = "tsne",split.by = 'group') 


###进行正常组和肿瘤组的差异分析
# 将group列赋值为Seurat对象的身份
##记住，后面分析要经常运行这行代码
yjsl <- JoinLayers(yjsl)
Idents(yjsl) <- yjsl$group

# 进行差异基因分析

# {
# DEGs <- FindMarkers(yjsl, ident.1 = 'PD', ident.2 = 'NC', logfc.threshold = 0.58)#可调节
# }
# 
# # 将分组信息添加到差异基因表格
# DEGs$group <- ifelse(DEGs$avg_log2FC > 0, "PD", "NC")  # 根据log2FoldChange决定分组
# 
# # 输出差异基因结果
# write.csv(DEGs, "PD和正常组差异基因.csv")
# 
# ## 过滤p值和调整后的p值
# DEGs.filter <- subset(DEGs, p_val_adj < 0.05)
# write.csv(DEGs.filter, "过滤后PD_vs_NC_差异基因.csv")

##重新赋予回来
Idents(yjsl) <- yjsl$seurat_clusters
###大体看一下情况###
# 查看多少个类型
table(yjsl@meta.data$seurat_clusters)
#只保留上调差异表达的基因
markers <- FindAllMarkers(yjsl, only.pos = TRUE, 
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
write.csv(markers,file = "cluster的标记基因.csv")
head(markers)
markers<-read.csv("cluster的标记基因.csv")
##保存用于后续分析
saveRDS(yjsl,"af2.rds")

library(randomcoloR)
##专门美化图片的包
#install.packages("/Users/zhangzhenhu/Desktop/单细胞注释_人工/scRNAtoolVis_0.1.0.tar.gz", repos=NULL, type="source")
# devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

yjsl=readRDS("af2.rds")

################推荐：人工手动注释################
######第一种，通过差异基因################

# 找到所有差异表达基因
# Idents(yjsl) <- yjsl$seurat_clusters
allmarkers <- markers#基因在至少 25% 的细胞中表达才会被考虑为差异基因
                            #只选择 log2 fold change 大于 0.25 的基因


# save(allmarkers,file = "allmarkers.rda")
# load("allmarkers.rda")

# # 提取每个簇中具有最显著差异表达的前4个基因
# top4_markers <- allmarkers %>%
#   group_by(cluster) %>%
#   top_n(n = 4, wt = avg_log2FC)
# 提取每个簇 avg_log2FC 最高的前100个基因，并输出CSV
top100_markers <- allmarkers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100) %>%
  ungroup() %>%
  dplyr::select(cluster, gene,
                avg_log2FC, pct.1, pct.2, p_val_adj)  # 加上 dplyr::

# 保存到CSV（用Excel打开）
write.csv(top100_markers, file = "top100_markers_per_cluster.csv", row.names = FALSE)

# #换个颜色的热图
# DoHeatmap(yjsl, features = unique(top4_markers$gene)) + NoLegend()+
#   scale_fill_gradientn(colors = c("#2fa1dd", "white", "#f87669"))
# #气泡图
# DotPlot(yjsl, features = unique(top4_markers$gene),cols = "RdYlBu") +
#   RotatedAxis()

###然后可以进一些注释网站
#CellMarker: http://xteam.xbio.top/CellMarker/

######第二种常见的标记基因################
# T Cells (CD3D, CD3E, CD8A)
# 
# B cells (CD19, CD79A, MS4A1 [CD20])
# 
# Plasma cells (IGHG1, MZB1, SDC1, CD79A)
# 
# Monocytes and macrophages (CD68, CD163, CD14)
# 
# NK Cells (FGFBP2, FCG3RA, CX3CR1)
# 
# Photoreceptor cells (RCVRN)
# 
# Fibroblasts (FGF7, MME)
# 
# Endothelial cells (PECAM1, VWF)
# 
# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24)
# 
# immune (CD45+,PTPRC)
# 
# epithelial/cancer (EpCAM+,EPCAM)
# 
# stromal (CD10+,MME,fibo or CD31+,PECAM1,endo)
# 
# mast cells( TPSAB1 and TPSB2 )
# 
# naive B cells(MS4A1 (CD20), CD19, CD22, TCL1A, and CD83)
# 
# plasma B cells(CD38, TNFRSF17 (BCMA), and IGHG1/IGHG4 


#######看别人的文献################


yjsl <- FindClusters(yjsl, resolution = 0.3)##记住和上面resolution数目一样

# genes <- list(
#   "T cell" = c("CD4", "CD3E", "CD8A"),
#   "B cell" = c("MS4A1", "CD79A", "PAX5"),
#   "Mono/Macro" = c("CD163", "CSF1R", "LYZ"),
#   "Plasma cell" = c("MZB1", "IGHG1", "XBP1"),
#   "Fibroblast" = c("COL12A1", "LUM", "COL1A1"),
#   "Endothelial" = c("PECAM1", "VWF", "AQP1"),
#   "Epithelial" = c("KRT6A", "KRT14", "DSG3"),
#   "pDC" = c("CLEC4C", "LILRA4", "GZMB"),
#   "MDSC" = c( "MKI67","TOP2A", "RRM2"),
#   "Neutrophils"=c("FCGR3B","CSF3R","G0S2"),
#   "mast cells"=c("TPSAB1","TPSB2","CPA3","KIT")
# )
#       
# do_DotPlot(sample = yjsl,features = genes,dot.scale = 10,legend.length = 50,
#            legend.framewidth = 2, font.size =10)
# 
# #另一种颜色气泡图
# DotPlot(yjsl, features = genes,cols = "RdYlBu") +
#   RotatedAxis()

#人工注释
table(yjsl$seurat_clusters)

table(yjsl@active.ident)
# 正确的细胞类型注释（按簇编号 0~14 顺序）
ann.ids <- c(
  "T cells",                  # 0 CD3D/CD3E/TRBC1/CCL5
  "Endothelial cells",        # 1 ACKR1/VWF/SELP/SELE inflammatory EC
  "Fibroblasts",              # 2 COL1A1/LUM/DCN/PDGFRA
  "Plasma cells",             # 3 IGHG/IGLC/IGKC
  "Macrophages",              # 4 C1QA/C1QB/CD163/CSF1R
  "Epithelial cells",         # 5 KRT14/KRT15/DSG3
  "Pericytes",                # 6 RGS5/ACTA2/TAGLN/MYH11
  "B cells",                  # 7 MS4A1/CD79A/PAX5
  "Endothelial cells",        # 8 GJA5/GJA4/SEMA3G arterial EC
  "Mast cells",               # 9 TPSAB1/CPA3/MS4A2
  "Neutrophils",              # 10 FCGR3B/CSF3R
  "Endothelial cells", # 11 CCL21/PROX1/PDPN/FLT4
  "Cycling cells",            # 12 MKI67/TOP2A
  "pDCs",                     # 14 CLEC4C/LILRA4
  "Pericytes",                # 15 RGS5/COX4I2/KCNJ8
  "Cycling cells",            # 16 BIRC5/PBK/UBE2C
  "Plasma cells"              # 17 MZB1/IGHG
)

afidens=mapvalues(Idents(yjsl), from = levels(Idents(yjsl)), to = ann.ids)
Idents(yjsl)=afidens
yjsl$cellType=Idents(yjsl)

#########人工注释后结果可视化
# 可视化UMAP/tSNE
# 查看当前簇编号（确认是13）
table(yjsl$seurat_clusters)
Idents(yjsl) <- yjsl$seurat_clusters
# # 删除簇13
# yjsl <- subset(yjsl, idents = 13, invert = TRUE)

# # 重新分配注释
# yjsl$celltype <- ann.ids[as.numeric(Idents(yjsl)) + 1]

DimPlot(yjsl, reduction = "umap", label = T, label.size = 4)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                                                                     legend.position = "right")

DimPlot(yjsl, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 3,
        split.by = "group") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        legend.position = "right")

# DimPlot(yjsl, reduction = "tsne", label = T, label.size = 3.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right",split.by = 'group')

saveRDS(yjsl, file = "yjsl_clean_allGSE.rds")
yjsl <- readRDS("yjsl_clean_allGSE.rds")

{
  Endothelial <- NormalizeData(Endothelial, normalization.method = "LogNormalize", scale.factor = 1e4) 
  #鉴定高变基
  Endothelial<- FindVariableFeatures(Endothelial, selection.method = 'vst', nfeatures = 2000)
  #数据Scaling，降维前的必要准备
  all.genes <- rownames(Endothelial)
  Endothelial <- ScaleData(Endothelial, features = all.genes)
  #正式降维
  Endothelial <- RunPCA(Endothelial, features = VariableFeatures(object = Endothelial)) 
}
#数据可视化方法
#VizDimReduction
VizDimLoadings(Endothelial, dims = 1:2, reduction = "pca")
#ElbowPlot
ElbowPlot(Endothelial, ndims=20, reduction="pca") 

########确定维度、细胞分群########
#选择PC
PC=1:15

###根据上一步中的合适维度对细胞进行分群
Endothelial <- FindNeighbors(Endothelial, dims = PC)
#resolution可以设0.1-1之间，值越高，亚群数目越多，常规0.5
Endothelial <- FindClusters(Endothelial, resolution =seq(0.3,0.8,0.1))
clustree(Endothelial)

Endothelial <- FindClusters(Endothelial, resolution =0.5)##记住这个维度
###使用UMAP及进行非线性降维
Endothelial <- RunUMAP(Endothelial, dims = PC)
#数据可视化
DimPlot(Endothelial, reduction = "umap", label = T,repel = T)

# ###使用tSNE进行非线性降维
# Endothelial = RunTSNE(Endothelial, dims = PC)
# DimPlot(Endothelial, reduction = "tsne") 

# 查看多少个类型
table(Endothelial@meta.data$seurat_clusters)
#只保留上调差异表达的基因
markers <- FindAllMarkers(Endothelial, only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)
write.csv(markers,file = "单独Endothelial的标记基因.csv")
head(markers)
# 提取每个簇 avg_log2FC 最高的前50个基因，并输出CSV
top100_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100) %>%
  ungroup() %>%
  dplyr::select(cluster, gene,
                avg_log2FC, pct.1, pct.2, p_val_adj)  # 加上 dplyr::

# 保存到CSV（用Excel打开）
write.csv(top100_markers, file = "top100_markers_EC_percluster.csv", row.names = FALSE)

#注释
{
  Endothelial <- FindClusters(Endothelial, resolution = 0.5)##记住和上面resolution数目一样
}

Endothelial_clean <- subset(
  Endothelial,
  idents = c(
    "0","5","6","7","9",
    "10","11","12","13","14"
  ),
  invert = TRUE
)
saveRDS(Endothelial_clean,"Endothelial_clean.rds")
# Endothelial <- subset(Endothelial, idents = "11", invert = TRUE)
Endothelial<-readRDS("Endothelial_clean.rds")
Endothelial <- NormalizeData(
  Endothelial,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

Endothelial <- FindVariableFeatures(
  Endothelial,
  selection.method = "vst",
  nfeatures = 2000
)

Endothelial <- ScaleData(Endothelial)

Endothelial <- RunPCA(
  Endothelial,
  features = VariableFeatures(Endothelial)
)

ElbowPlot(Endothelial)

Endothelial <- FindNeighbors(
  Endothelial,
  dims = 1:13
)

Endothelial <- FindClusters(
  Endothelial,
  resolution = 0.5
)

Endothelial <- RunUMAP(
  Endothelial,
  dims = 1:13
)

DimPlot(
  Endothelial,
  reduction="umap",
  label=TRUE
)

markers <- FindAllMarkers(Endothelial, only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)
write.csv(markers,file = "单独Endothelial的标记基因2.csv")
Idents(Endothelial) <- Endothelial$seurat_clusters

Endothelial_clean <- subset(
  Endothelial,
  idents = "5",
  invert = TRUE
)
ann.ids <- c(
  "VEC",
  "CVEC",
  "Inflammatory_EC",
  "AVEC",
  "Activated_EC",
  "NA",
  "LEC",
  "LEC"
)

# 4. 将注释信息赋给 Idents 和 meta.data
Idents(Endothelial_clean) <- Endothelial_clean$seurat_clusters
new_ids <- ann.ids[as.numeric(as.character(Idents(Endothelial_clean))) + 1]  # 因为 cluster 编号从0开始
Endothelial_clean$subtype <- new_ids
Idents(Endothelial_clean) <- new_ids


table(Endothelial_clean$subtype)
table(Endothelial_clean$seurat_clusters)



DimPlot(
  Endothelial_clean,
  reduction="umap",
  group.by="subtype",
  label=TRUE,
  repel=TRUE
)+
  theme_classic()

colaa=distinctColorPalette(100)

library(ggplot2)
library(Seurat)
library(dplyr)
library(ggrepel)


# ===============================
# color palette
# ===============================

ec_colors <- c(
  
  "Activated_EC"="#E64B35",
  "Inflammatory_EC"="#4DBBD5",
  "AVEC"="#00A087",
  "CVEC"="#3C5488",
  "LEC"="#8491B4",
  "VEC"="#F39B7F"
  
)



group_colors <- c(
  "NC"="#999999",
  "PD"="#E64B35"
)


theme_sci <- theme_classic()+
  theme(
    
    plot.title=
      element_text(
        size=16,
        face="bold",
        hjust=0.5
      ),
    
    axis.title=
      element_text(
        size=14
      ),
    
    axis.text=
      element_text(
        size=12
      ),
    
    legend.title=
      element_text(
        size=13,
        face="bold"
      ),
    
    legend.text=
      element_text(
        size=11
      )
    
  )

Endothelial_clean$subtype <- factor(
  Endothelial_clean$subtype,
  levels=c(
    "Activated_EC",
    "Inflammatory_EC",
    "AVEC",
    "CVEC",
    "LEC",
    "VEC"
  )
)


p_marker <- DotPlot(
  Endothelial_clean,
  features=endo_markers2,
  group.by="subtype",
  cols=c(
    "lightgrey",
    "#3C5488"
  )
)+
  RotatedAxis()+
  theme_sci+
  ggtitle(
    "Marker expression defines endothelial subtypes"
  )

p_marker

VlnPlot(
  Endothelial_clean,
  features = "Activated_signature1",
  group.by = "subtype",
  cols = ec_colors,
  pt.size = 0
) +
  theme_sci +
  ggtitle("Activated endothelial transcriptional program") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



Endothelial_clean<-read_rds("Endothelial_clean_final.rds")
activated_program <- c(
  "ARL4C",
  "NFKBIZ",
  "SERPINE1",
  "C2CD4B",
  "BHLHE40",
  "SELP",
  "SELE",
  "ICAM1"
)

inflam_markers <- FindMarkers(
  Endothelial_clean,
  ident.1="Inflammatory_EC",
  min.pct=0.25,
  logfc.threshold=0.25
)


write.csv(
  inflam_markers,
  "Inflammatory_EC_DEG.csv"
)

head(inflam_markers,30)

inflammatory_program <- c(
  "ACKR1",
  "IL1R1",
  "HLA-DPA1",
  "HLA-DRA",
  "HLA-DRB5",
  "CSF2RB",
  "VCAM1",
  "SELE"
)
Endothelial_clean <- AddModuleScore(
  Endothelial_clean,
  features=list(activated_program),
  name="Activated_program"
)


Endothelial_clean <- AddModuleScore(
  Endothelial_clean,
  features=list(inflammatory_program),
  name="Inflammatory_program"
)

VlnPlot(
  Endothelial_clean,
  features=c(
    "Activated_program1",
    "Inflammatory_program1"
  ),
  group.by="subtype",
  pt.size=0,
  combine=FALSE
)

library(ggplot2)
library(ggpubr)

score_df <- FetchData(
  Endothelial_clean,
  vars=c(
    "Activated_score1",
    "Inflammatory_score1",
    "subtype"
  )
)


cor_test <- cor.test(
  score_df$Activated_score1,
  score_df$Inflammatory_score1,
  method="spearman"
)


cor_test

library(ComplexHeatmap)
library(circlize)


avg_score <- score_df %>%
  dplyr::group_by(subtype) %>%
  dplyr::summarise(
    Activated_program =
      mean(Activated_program1),
    Inflammatory_program =
      mean(Inflammatory_program1)
  )


heat_data <- avg_score %>%
  tibble::column_to_rownames("subtype")


heat_data <- as.matrix(heat_data)


Heatmap(
  heat_data,
  name="Module score",
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  column_title="Endothelial activation programs",
  row_names_gp=gpar(fontsize=12),
  column_names_gp=gpar(fontsize=12),
  col=colorRamp2(
    c(min(heat_data),
      0,
      max(heat_data)),
    c("#4575B4",
      "white",
      "#D73027")
  )
)
#####以后的代码####

#####CSF2####
#最常规的小提琴图
VlnPlot(yjsl, features = "CSF2",group.by="cellType")+NoLegend()

table(yjsl$group)

library(ggplot2)
DotPlot(yjsl, features = "CSF2", group.by = "cellType") +
  coord_flip() +   # 横向显示，方便阅读
  labs(title = "CSF2 expression across cell types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(data.table)
library(ggplot2)

# 1. 提取数据（确保为 data.table）
expr_data <- FetchData(yjsl, vars = c("CSF2", "cellType", "group")) %>%
  na.omit() %>%
  as.data.table()

# 2. 按 cellType 和 group 分组计算统计量
summ <- expr_data[, .(
  avg_exp = mean(CSF2, na.rm = TRUE),
  pct_pos = mean(CSF2 > 0, na.rm = TRUE) * 100
), by = .(cellType, group)]

# 检查结果
print(summ)

# 3. 按平均表达排序细胞类型
mean_per_cell <- summ[, .(mean_overall = mean(avg_exp)), by = cellType]
overall_order <- mean_per_cell[order(-mean_overall)]$cellType
summ$cellType <- factor(summ$cellType, levels = overall_order)

# 4. 绘图
ggplot(summ, aes(x = cellType, y = group, size = pct_pos, color = avg_exp)) +
  geom_point() +
  scale_size_continuous(range = c(0, 10), name = "% positive") +
  scale_color_gradient(low = "#00BFC4", high = "salmon", name = "Avg. exp") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Cell Type", y = "", title = "CSF2 expression (PD vs NC)")

# 提取数据
plot_data <- FetchData(yjsl, vars = c("CSF2", "cellType"))
plot_data <- plot_data[plot_data$CSF2 > 0, ]   # 只保留表达>0的细胞

ggplot(plot_data, aes(x = cellType, y = CSF2)) +
  geom_jitter(aes(color = cellType), width = 0.2, size = 0.6, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CSF2-expressing cells only",
       y = "Log-normalized expression", x = "")


p1<-FeaturePlot(yjsl,features = "CSF2")
p1


FeaturePlot(yjsl, features = "CSF2", max.cutoff = "q75", pt.size = 0.3 ,order = TRUE) +
  ggtitle("CSF2 expression")


library(data.table)
library(ggplot2)

# 1. 提取数据（确保为 data.table）
expr_data <- FetchData(yjsl_filtered, vars = c("CYP19A1", "cellType", "group")) %>%
  na.omit() %>%
  as.data.table()

# 2. 按 cellType 和 group 分组计算统计量
summ <- expr_data[, .(
  avg_exp = mean(CYP19A1, na.rm = TRUE),
  pct_pos = mean(CYP19A1 > 0, na.rm = TRUE) * 100
), by = .(cellType, group)]

# 检查结果
print(summ)

# 3. 按平均表达排序细胞类型
mean_per_cell <- summ[, .(mean_overall = mean(avg_exp)), by = cellType]
overall_order <- mean_per_cell[order(-mean_overall)]$cellType
summ$cellType <- factor(summ$cellType, levels = overall_order)

# 4. 绘图
ggplot(summ, aes(x = cellType, y = group, size = pct_pos, color = avg_exp)) +
  geom_point() +
  scale_size_continuous(range = c(0, 10), name = "% positive") +
  scale_color_gradient(low = "#00BFC4", high = "salmon", name = "Avg. exp") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Cell Type", y = "", title = "CYP19A1 expression (PD vs NC)")


# 查看每个类型有多少细胞
table(yjsl$cellType)
#####表达量分析#####
# 加载必要的包
# 确保 expr_data 已经准备好
expr_data <- FetchData(yjsl, vars = c("CSF2", "group", "cellType"))
expr_data <- expr_data %>% filter(!is.na(group))

cell_types <- unique(expr_data$cellType)
results <- data.frame()

for (ct in cell_types) {
  sub_data <- expr_data %>% filter(cellType == ct)
  
  if (sum(sub_data$group == "NC") > 1 & sum(sub_data$group == "PD") > 1) {
    # 注意：此处将 COX4I1 改为 CSF2
    mean_NC <- mean(sub_data$CSF2[sub_data$group == "NC"])
    mean_PD <- mean(sub_data$CSF2[sub_data$group == "PD"])
    log2FC <- log2((mean_PD + 0.01) / (mean_NC + 0.01))
    test_result <- wilcox.test(CSF2 ~ group, data = sub_data)
    temp_row <- data.frame(
      cellType = ct,
      mean_NC = mean_NC,
      mean_PD = mean_PD,
      log2FC = log2FC,
      p_value = test_result$p.value,
      n_NC = sum(sub_data$group == "NC"),
      n_PD = sum(sub_data$group == "PD")
    )
    results <- rbind(results, temp_row)
  } else {
    temp_row <- data.frame(
      cellType = ct,
      mean_NC = NA,
      mean_PD = NA,
      log2FC = NA,
      p_value = NA,
      n_NC = sum(sub_data$group == "NC"),
      n_PD = sum(sub_data$group == "PD")
    )
    results <- rbind(results, temp_row)
  }
}



# ============================================================================
# 细胞类型标记基因表达点图（DotPlot）- 阶梯状排列
# 基于 yjsl 对象（已有人工注释的 cellType）
# ============================================================================
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# 1. 定义细胞类型及其标记基因（与您的 genes 列表一致，并调整名称匹配 cellType）
# 注意：将 "Neutrophils" 改为 "Neutrophil" 以匹配您的注释
# 同样 "mast cells" 改为 "Mast cell"
lineage_markers<- list(
  
  # T cells
  "T cells" = c(
    "CD3D",
    "CD3E"
    # "TRBC1",
    # "CD8A"
  ),
  
  # B cells
  "B cells" = c(
    "MS4A1",
    # "CD79A"
    "CD37"
    # "CD74"
  ),
  
  # Plasma cells
  "Plasma cells" = c(
    "MZB1",
    # "JCHAIN",
    "IGHG1"
    # "SDC1"
  ),
  
  # Macrophages
  "Macrophages" = c(
       "LYZ","C1QA"
    # "C1QB",
    # "CD163",
 
  ),
  
  # Fibroblasts
  "Fibroblasts" = c(
    "COL1A1",
    "COL1A2"
    # "DCN",
    # "LUM"
  ),
  
  # Endothelial cells
  "Endothelial cells" = c(
    "PECAM1",
    "VWF","CDH5"
    # "EMCN"
  ),
  
  # Pericytes
  "Pericytes" = c(
    "RGS5",
    "ACTA2"
    # "TAGLN",
    # "MCAM"
  ),
  
  # Epithelial cells
  "Epithelial cells" = c(
    "KRT14",
    "KRT5"
    # "KRT19",
    # "EPCAM"
  ),
  
  # Mast cells
  "Mast cells" = c(
    "TPSAB1",
    "TPSB2"
    # "CPA3",
    # "KIT"
  ),
  
  # Neutrophils
  "Neutrophils" = c(
    "FCGR3B",
    "CSF3R"
    # "S100A8",
    # "S100A9"
  ),
  
  # pDC
  "pDCs" = c(
    # "CLEC4C",
    "GZMB",
    "IRF7"
    # "IL3RA"
  ),
  
  # Cycling
  "Cycling cells" = c(
    "MKI67",
    "TOP2A"
    # "BIRC5"
  )
  
)

# ==================== 过滤只在数据中存在的基因 ====================
available_markers <- list()
for(cell_type in names(lineage_markers)) {
  genes <- lineage_markers[[cell_type]]
  present <- genes[genes %in% rownames(yjsl)]
  if(length(present) > 0) {
    available_markers[[cell_type]] <- present
  } else {
    message("警告：细胞类型 ", cell_type, " 的所有标记基因均未在数据中找到！")
  }
}

# 打印实际使用的基因
cat("=== 实际使用的标记基因 ===\n")
for(cell_type in names(available_markers)) {
  cat(cell_type, ": ", paste(available_markers[[cell_type]], collapse = ", "), "\n")
}

# ==================== 确保 cellType 列存在且按顺序排列 ====================
if (!"cellType" %in% colnames(yjsl@meta.data)) {
  stop("yjsl 中缺少 cellType 列，请先运行注释代码！")
}

# 只保留 available_markers 中出现的细胞类型（过滤掉不在列表中的类型，避免报错）
valid_types <- intersect(names(available_markers), unique(yjsl$cellType))
if(length(valid_types) == 0) {
  stop("没有匹配的细胞类型！请检查 cellType 列内容是否与列表名称一致。")
}

# 将 cellType 设为因子，顺序按列表（确保纵坐标顺序）
yjsl$cellType <- factor(yjsl$cellType, levels = valid_types)

# ==================== 绘制 DotPlot ====================
# 仅保留 available_markers 中对应的类型
features_plot <- available_markers[valid_types]

p <- DotPlot(yjsl, 
             features = features_plot,
             group.by = "cellType",
             cols = c("lightgrey", "steelblue"),
             dot.scale = 8) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
   
    x = "Marker genes",
    y = "Cell types"
  ) +
  scale_color_gradient2(
    low = "#00BFC4", 
    mid = "white", 
    high = "#F8766D",
    midpoint = 1,
    name = "Average Expression"
  )

# 显示图形
print(p)

# ==================== 保存图片 ====================
ggsave("DotPlot_celltype_markers.pdf", p, width = 16, height = 4)
