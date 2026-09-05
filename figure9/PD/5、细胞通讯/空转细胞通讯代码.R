library(Seurat)
library(CellChat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

## 读取空间转录组对象
load("牙周炎_注释.rda")


##  设置分组信息
group_col <- "Region"## group_col 表示用于 CellChat 分组的列

## 查看每个 Region 中有多少 spot
table(stRNA@meta.data[[group_col]])
  
# 定义实际的分组名称（与 Region 列一致）
group_levels <- c("Suprabasal keratinized",   # 表层角化
                  "Suprabasal granular",      # 颗粒层
                  "Suprabasal spinous",       # 棘层
                  "Basal epithelium",         # 基底层
                  "Junctional epithelium",    # 结合上皮（牙周炎特有）
                  "Immune stroma",            # 免疫基质（浅层基质）
                  "Reticular stroma"          # 网状基质（深层基质）
                  )

# 为每个组分配颜色（可以沿用您之前注释时用的颜色）
group_colors <- c(
  "Suprabasal keratinized"  = "#F7D9B0",  # 浅米杏（表层角化）
  "Suprabasal granular"     = "#EAA369",  # 暖杏橙（颗粒层）
  "Suprabasal spinous"      = "#D99B5A",  # 深杏/橘棕（棘层）
  "Basal epithelium"        = "#EF9B94",  # 鲑鱼粉（基底层，与健康一致）
  "Junctional epithelium"   = "#D35D5D",  # 砖红（结合上皮，炎症特征）
  "Immune stroma"           = "#6B8DAE",  # 冷灰蓝（免疫基质，与健康一致）
  "Reticular stroma"        = "#72B3AD"   # 冷蓝绿（网状基质，与健康一致）
)

##  提取表达矩阵
## CellChat 需要输入标准化后的表达矩阵，一般使用 Seurat assay 中的 data 层
{DefaultAssay(stRNA) <- "Spatial"
  
  data.input <- GetAssayData(
    stRNA,
    assay = "Spatial",
    layer = "data"
  )
}

## 如果 data 层为空，说明还没有 NormalizeData
## 这时先进行标准化
if (nrow(data.input) == 0 || ncol(data.input) == 0) {
  
  stRNA <- NormalizeData(stRNA)
  
  data.input <- GetAssayData(
    stRNA,
    assay = "Spatial",
    layer = "data"
  )
}

##  构建 CellChat 所需 meta 信息
{meta <- data.frame(
  labels = stRNA@meta.data[[group_col]],
  row.names = colnames(stRNA)
)
  
  ## 去掉没有 Region 注释的 spot
  meta <- meta[!is.na(meta$labels), , drop = FALSE]
}

## 设置 Region 顺序
meta$labels <- factor(meta$labels,levels = group_levels)

## 表达矩阵也只保留有 Region 注释的 spot
data.input <- data.input[, rownames(meta)]

## 创建 CellChat 对象
{cellchat <- createCellChat(
  object = data.input,## object = 表达矩阵
  meta = meta,## meta = 分组信息
  group.by = "labels"## group.by = meta 中用于分组的列名
)
  
  ## 固定分组顺序，避免后面画图顺序乱
  cellchat@idents <- factor(
    cellchat@idents,
    levels = group_levels
  )
}

## 设置配体-受体数据库， 使用人类 CellChat 数据库
CellChatDB <- CellChatDB.human

## 只使用分泌型信号通路，例如：PTN、SPP1、MIF、CXCL、CCL 等
CellChatDB.use <- subsetDB(
  CellChatDB,
  search = "Secreted Signaling"
)

## 将数据库加入 CellChat 对象
{cellchat@DB <- CellChatDB.use
  
  ## CellChat 标准分析流程
  ## subsetData： 从表达矩阵中提取数据库相关的配体和受体基因
  cellchat <- subsetData(cellchat)
  
  ## identifyOverExpressedGenes：找到每个 Region 中高表达的配体/受体基因
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  ## identifyOverExpressedInteractions：根据高表达基因筛选可能存在的配体-受体互作
  cellchat <- identifyOverExpressedInteractions(cellchat)
}

## computeCommunProb：计算通讯概率
{cellchat <- computeCommunProb(
  cellchat,
  type = "truncatedMean",## type = "truncatedMean" 表示使用截尾均值，减少极端值影响
  trim = 0.1## trim = 0.1 表示去掉两端10%的极端表达值
)
  
  
  ## filterCommunication：过滤细胞数太少的分组
  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )
}

#保存
save(cellchat,file = "cellchat.rda")
load("cellchat.rda")

## computeCommunProbPathway：将单个配体-受体对整合到信号通路层面
cellchat <- computeCommunProbPathway(cellchat)

## aggregateNet：汇总所有通讯网络，得到通讯数量 count 和通讯强度 weight
cellchat <- aggregateNet(cellchat)

## 计算每个 Region 在通讯网络中的中心性，用于判断谁是主要发送者、接收者、中介者等
{cellchat <- netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)
}

## 提取所有通讯结果，subsetCommunication 会导出所有显著配体-受体通讯结果
df.net <- subsetCommunication(cellchat)

write.csv(
  df.net,
  "CellChat_all_communication_results.csv",
  row.names = FALSE
)


## 通讯数量圈图， groupSize 表示每个 Region 的 spot 数量
groupSize <- as.numeric(table(cellchat@idents))
par(
  mar = c(2,2,4,2),
  xpd = TRUE
)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 1.1,
  vertex.size.max = 25,
  color.use = group_colors,
  title.name = "Number of Interactions"
)

## 通讯强度圈图
netVisual_circle(
  cellchat@net$weight,      # 通讯强度矩阵
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 1.1,
  vertex.size.max = 25,
  color.use = group_colors,
  title.name = "Interaction Weights/Strength"
)

## 通讯数量热图
netVisual_heatmap(
  cellchat,
  measure = "count",
  color.heatmap = "Blues",
  font.size = 12,
  font.size.title = 14,
  title.name = "Number of Interactions"
)

##通讯强度热图

netVisual_heatmap(
  cellchat,
  measure = "weight",
  color.heatmap = "Blues",
  font.size = 12,
  font.size.title = 14,
  title.name = "Interaction Weights/Strength"
)

## 查看检测到的信号通路
pathways.all <- cellchat@netP$pathways

print(pathways.all)

write.csv(
  data.frame(pathway = pathways.all),
  "CellChat_detected_pathways.csv",
  row.names = FALSE
)
# # 提取通路强度
# pathway_strength_spatial <- apply(cellchat@netP$prob, 3, sum)
# sorted_spatial <- sort(pathway_strength_spatial, decreasing = TRUE)
# print(sorted_spatial)

# 提取所有通路的发送者-接收者强度
netP <- cellchat@netP$prob
# 计算每个通路在每组的总发送强度
sender_matrix <- apply(netP, c(1, 3), sum)   # Region × Pathway
# 绘制热图
pheatmap(sender_matrix,
         main = "Sender strength per region",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "steelblue", "red"))(100))



##  筛选重点关注通路
## 这里列出一些肿瘤微环境中常见通路
selected_pathways <- intersect(
  c(
    "MIF",
    "SELE",
    "PECAM1",
    "CDH5",
    "COLLAGEN",
    "MIF",
    "CXCL"
  ),
  pathways.all
)

print(selected_pathways)

## 信号通路整体强度排序
## 回答：在region 之间，哪些信号通路整体最强？
df.pathway <- data.frame(
  pathway = cellchat@netP$pathways,
  total_prob = apply(cellchat@netP$prob, 3, sum)
) %>%
  arrange(desc(total_prob))

write.csv(
  df.pathway,
  "Pathway_total_strength_rank.csv",
  row.names = FALSE
)

p_pathway_rank <- ggplot(
  df.pathway[1:20, ],
  aes(
    x = reorder(pathway, total_prob),
    y = total_prob
  )
) +
  geom_col(fill = "#3B5BA5") +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Total communication probability",
    title = "Top signaling pathways"
  )

p_pathway_rank
# PERIOSTIN
# TGFb
# PDGF
########PERIOSTIN####
pathway.show <- "PERIOSTIN"

par(
  mar = c(1, 1, 3, 1),
  xpd = TRUE
)
{
  netVisual_aggregate(
    cellchat,
    signaling = pathway.show,
    layout = "circle",
    color.use = group_colors,
    vertex.label.cex = 1.1,
    vertex.size.max = 20,
    edge.width.max = 10,
    arrow.width = 1.2,
    arrow.size = 0.5,
    top = 0.8
  )
}


## 提取cxcl 通路通讯矩阵
net.cxcl <- cellchat@netP$prob[, , pathway.show]

net.cxcl[is.na(net.cxcl)] <- 0
net.cxcl

# 提取原始坐标
coords <- GetTissueCoordinates(stRNA)
if (all(c("x", "y") %in% colnames(coords))) {
  coords_use <- coords[, c("x", "y")]
} else if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
  coords_use <- coords[, c("imagecol", "imagerow")]
  colnames(coords_use) <- c("x", "y")
} else {
  stop("坐标列名不匹配")
}

# 保留行名
rownames_orig <- rownames(coords_use)

# 变换：顺时针90° + 左右翻转 → 等效于 x' = -y, y' = -x
coords_use <- data.frame(
  x = coords_use$y,
  y = coords_use$x
)
rownames(coords_use) <- rownames_orig

# 现在按 stRNA 的 spot 顺序排列（这会保留行名）
coords_use <- coords_use[colnames(stRNA), ]

# 后续步骤不变...
meta_spatial <- data.frame(
  labels = stRNA$Region,
  row.names = colnames(stRNA)
)
meta_spatial <- meta_spatial[rownames(coords_use), , drop = FALSE]
meta_spatial$labels <- factor(meta_spatial$labels, 
                              levels = c("Suprabasal keratinized",
                                         "Suprabasal granular",
                                         "Suprabasal spinous",
                                         "Basal epithelium",
                                         "Junctional epithelium",
                                         "Immune stroma",
                                         "Reticular stroma"))


# 4. 检查各区域数量
table(meta_spatial$labels)

## cxcl空间通讯图
par(
  mar = c(1, 1, 3, 1),
  xpd = TRUE
)

netVisual_spatial(
  net = net.cxcl,
  coordinates = coords_use,
  labels = meta_spatial$labels,
  color.use = group_colors,
  vertex.size.max = 3,        # 调大顶点
  vertex.weight = 40,          # 调大权重
  vertex.label.cex = 4,
  edge.width.max = 20,
  alpha.edge = 1,
  point.size = 3,
  alpha.image = 0.7,
  remove.isolate = FALSE,
  title.name = "PERIOSTIN signaling spatial network"
)



## PERIOSTIN Bubble 图， Bubble 图可以显示：哪些 Region 对之间存在 PTN 通讯
## 点大小/颜色通常代表通讯强度或显著性

netVisual_bubble(
  cellchat,
  sources.use = c("Immune stroma"),
  targets.use = c("Suprabasal keratinized",
                  "Suprabasal granular",
                  "Suprabasal spinous",
                  "Basal epithelium",
                  "Junctional epithelium",
                  "Immune stroma",
                  "Reticular stroma"),
  signaling = pathway.show,
  remove.isolate = FALSE,
  angle.x = 90,
  color.text = FALSE,
  font.size = 12,
  line.on = TRUE
  
)

## CXCL chord gene 图
## chord图展示：哪些配体-受体对参与 PTN 通路通讯

circos.clear()
par(
  mar = c(1, 1, 3, 1)
)

netVisual_chord_gene(
  cellchat,
  signaling = pathway.show,
  lab.cex = 1.2,
  big.gap = 20,
  small.gap = 3,
  reduce = 0.05,
  transparency = 0.15,
  color.use = group_colors
)

title(
  main = "PERIOSTIN signaling network",
  cex.main = 1.5,
  font.main = 2
)

## cxcl 通路贡献图，显示cxcl通路中哪些配体-受体对贡献最大
print(
  netAnalysis_contribution(
    cellchat,
    signaling = pathway.show
  )
)

## SPP1 通路相关基因在 Region 中的表达小提琴图
genes.use <- c(
  "CXCL12",
  "CXCL13",
  "CXCL2",
  "CXCL6",
  "ACKR3",
  "CXCR4",
  "ACKR1",
  "CXCR2"
  )

genes.use <- intersect(
  genes.use,
  rownames(stRNA)
)

p_gene_vln <- VlnPlot(
  stRNA,
  features = genes.use,
  group.by = "Region",
  pt.size = 0,
  ncol = 4
)

p_gene_vln


########TGFb####
## TGFb通路网络图
## 先确认 mif是否存在
pathway.show <- "TGFb"
par(
  mar = c(1, 1, 3, 1),
  xpd = TRUE
)
{
  netVisual_aggregate(
    cellchat,
    signaling = pathway.show,
    layout = "circle",
    color.use = group_colors,
    vertex.label.cex = 1.1,
    vertex.size.max = 10,
    edge.width.max = 10,
    arrow.width = 1.2,
    arrow.size = 0.5,
    top = 0.8
  )
}

# 提取原始坐标
coords <- GetTissueCoordinates(stRNA)
if (all(c("x", "y") %in% colnames(coords))) {
  coords_use <- coords[, c("x", "y")]
} else if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
  coords_use <- coords[, c("imagecol", "imagerow")]
  colnames(coords_use) <- c("x", "y")
} else {
  stop("坐标列名不匹配")
}

# 保留行名
rownames_orig <- rownames(coords_use)

# 变换：顺时针90° + 左右翻转 → 等效于 x' = -y, y' = -x
coords_use <- data.frame(
  x = coords_use$y,
  y = coords_use$x
)
rownames(coords_use) <- rownames_orig

# 现在按 stRNA 的 spot 顺序排列（这会保留行名）
coords_use <- coords_use[colnames(stRNA), ]

##  整理空间 meta 信息
meta_spatial <- data.frame(
  labels = stRNA$Region,
  row.names = colnames(stRNA)
)

## 保证 meta 和 coords 顺序一致
meta_spatial <- meta_spatial[rownames(coords_use), , drop = FALSE]

meta_spatial$labels <- factor(
  meta_spatial$labels,
  levels = c("Suprabasal keratinized",   # 表层角化
             "Suprabasal granular",      # 颗粒层
             "Suprabasal spinous",       # 棘层
             "Basal epithelium",         # 基底层
             "Junctional epithelium",    # 结合上皮（牙周炎特有）
             "Immune stroma",            # 免疫基质（浅层基质）
             "Reticular stroma"          # 网状基质（深层基质）
  )
)

all(rownames(coords_use) == rownames(meta_spatial))

## 提取MIF 通路通讯矩阵
net.MIF <- cellchat@netP$prob[, , pathway.show]

net.MIF[is.na(net.MIF)] <- 0
net.MIF

# 4. 检查各区域数量
table(meta_spatial$labels)

## MIF空间通讯图
par(
  mar = c(1, 1, 3, 1),
  xpd = TRUE
)

netVisual_spatial(
  net = net.MIF,
  coordinates = coords_use,
  labels = meta_spatial$labels,
  color.use = group_colors,
  vertex.size.max = 4,        # 调大顶点
  vertex.weight = 60,          # 调大权重
  vertex.label.cex = 4,
  edge.width.max = 30,
  alpha.edge = 1,
  point.size = 3,
  alpha.image = 0.7,
  remove.isolate = FALSE,
  title.name = "TGFb signaling spatial network"
)

## MIF Bubble 图， Bubble 图可以显示：哪些 Region 对之间存在 PTN 通讯
## 点大小/颜色通常代表通讯强度或显著性

netVisual_bubble(
  cellchat,
  sources.use = c("Immune stroma"),
  targets.use = c("Suprabasal keratinized",
                  "Suprabasal granular",
                  "Suprabasal spinous",
                  "Basal epithelium",
                  "Junctional epithelium",
                  "Immune stroma",
                  "Reticular stroma"),
  signaling = pathway.show,
  remove.isolate = FALSE,
  angle.x = 90,
  color.text = FALSE,
  font.size = 12,
  line.on = TRUE
  
)

## MIF chord gene 图
## chord图展示：哪些配体-受体对参与 PTN 通路通讯

circos.clear()
par(
  mar = c(1, 1, 3, 1)
)

netVisual_chord_gene(
  cellchat,
  signaling = pathway.show,
  lab.cex = 1.2,
  big.gap = 20,
  small.gap = 3,
  reduce = 0.05,
  transparency = 0.15,
  color.use = group_colors
)

title(
  main = "MIF signaling network",
  cex.main = 1.5,
  font.main = 2
)

## MIF通路贡献图，显示MIF通路中哪些配体-受体对贡献最大
print(
  netAnalysis_contribution(
    cellchat,
    signaling = pathway.show
  )
)

## MIF通路相关基因在 Region 中的表达小提琴图
genes.use <- c(
  "MIF",
  "CD44",
  "CD74",
  "ACKR3",
  "CXCR4",
  "CXCR2"
)

genes.use <- intersect(
  genes.use,
  rownames(stRNA)
)

p_gene_vln <- VlnPlot(
  stRNA,
  features = genes.use,
  group.by = "Region",
  pt.size = 0,
  ncol = 3
)

p_gene_vln
