############################################################
# AGT virtual overexpression in PD Fibroblasts
############################################################

library(scTenifoldNet)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(Matrix)

set.seed(12345)

############################################################
# 1. 工作目录和数据读取
############################################################

setwd("G:/1Yunvjian/0a26.7.7singlecell/oe/csf2")

yjsl <- readRDS("yjsl_clean_allGSE.rds")

DefaultAssay(yjsl) <- "RNA"

############################################################
# 2. 检查元数据字段和细胞类型
############################################################

required_meta <- c("cellType", "group")

missing_meta <- setdiff(
  required_meta,
  colnames(yjsl@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "缺少以下元数据字段：",
    paste(missing_meta, collapse = ", ")
  )
}

if (!"Fibroblasts" %in% unique(as.character(yjsl$cellType))) {
  cat("当前cellType包括：\n")
  print(sort(unique(as.character(yjsl$cellType))))
  stop("cellType中不存在Fibroblasts，请根据上面打印的名称检查注释")
}

if (!"PD" %in% unique(as.character(yjsl$group))) {
  cat("当前group包括：\n")
  print(sort(unique(as.character(yjsl$group))))
  stop("group中不存在PD")
}

############################################################
# 3. 提取PD组Fibroblasts
############################################################

oe_obj <- subset(
  yjsl,
  subset = cellType == "Fibroblasts" & group == "PD"
)

cat("PD Fibroblasts cell number:", ncol(oe_obj), "\n")
print(dim(oe_obj))

if (ncol(oe_obj) < 200) {
  stop("PD Fibroblasts细胞数少于200，不建议进行scTenifoldOE")
}

############################################################
# 4. 提取RNA counts矩阵
############################################################

countMat <- GetAssayData(
  oe_obj,
  assay = "RNA",
  layer = "counts"
)

cat("Original count matrix dimension:\n")
print(dim(countMat))

############################################################
# 5. 检查AGT
############################################################

target_gene <- "AGT"

if (!target_gene %in% rownames(countMat)) {
  stop("AGT不存在于RNA counts矩阵")
}

AGT_cells <- sum(countMat[target_gene, ] > 0)

AGT_percent <- mean(countMat[target_gene, ] > 0)

cat(
  "AGT positive cells:",
  AGT_cells,
  "\n"
)

cat(
  "AGT expression fraction:",
  round(AGT_percent, 3),
  "\n"
)

############################################################
# AGT强制保留版本
############################################################

target_gene <- "AGT"

############################################################
# 1. 提取counts矩阵
############################################################

countMat <- GetAssayData(
  oe_obj,
  assay = "RNA",
  layer = "counts"
)

# 确保基因名存在
if (is.null(rownames(countMat))) {
  stop("countMat没有基因行名")
}

# 严格检查AGT是否存在于原始counts矩阵
if (!target_gene %in% rownames(countMat)) {
  cat("与AGT名称相似的基因：\n")
  print(
    grep(
      "AGT",
      rownames(countMat),
      value = TRUE,
      ignore.case = TRUE
    )
  )
  
  stop(
    "AGT不存在于原始RNA counts矩阵中，不能直接强制保留；",
    "请检查矩阵使用的是Gene Symbol还是Ensembl ID"
  )
}

# 使用drop=FALSE，防止单行矩阵降维和丢失行名
AGT_row <- countMat[
  target_gene,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(AGT_row),
    target_gene
  )
)

cat("AGT total counts:", sum(AGT_row), "\n")
cat("AGT positive cells:", sum(AGT_row > 0), "\n")
cat(
  "AGT expression fraction:",
  round(mean(AGT_row > 0), 4),
  "\n"
)

if (sum(AGT_row) == 0) {
  warning(
    "AGT虽然存在于矩阵中，但在当前PD Fibroblasts中全部为0；",
    "这种情况下网络无法可靠推断AGT的调控关系"
  )
}

############################################################
# 2. 基因表达过滤：明确强制保留AGT
############################################################

expr_ratio <- Matrix::rowSums(countMat > 0) / ncol(countMat)

keep_by_expression <- expr_ratio >= 0.02

# 不论AGT表达比例是否达到2%，都保留AGT
keep_by_expression[
  rownames(countMat) == target_gene
] <- TRUE

count_filtered <- countMat[
  keep_by_expression,
  ,
  drop = FALSE
]

# 再次检查
if (!target_gene %in% rownames(count_filtered)) {
  stop("AGT在表达过滤后丢失")
}

cat(
  "AGT retained after expression filtering:",
  target_gene %in% rownames(count_filtered),
  "\n"
)

############################################################
# 3. 选择高变基因
############################################################

oe_obj <- FindVariableFeatures(
  oe_obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 5000,
  verbose = FALSE
)

hvgs <- VariableFeatures(oe_obj)

# AGT必须放在基因列表首位
keep_gene <- unique(
  c(
    target_gene,
    hvgs
  )
)

keep_gene <- keep_gene[
  keep_gene %in% rownames(count_filtered)
]

# 使用match和drop=FALSE，确保行名及顺序不丢失
input_index <- match(
  keep_gene,
  rownames(count_filtered)
)

input_index <- input_index[
  !is.na(input_index)
]

data <- count_filtered[
  input_index,
  ,
  drop = FALSE
]

############################################################
# 4. 最终强制检查AGT
############################################################

if (!target_gene %in% rownames(data)) {
  
  # 理论上不会进入这里；作为额外保险再次加入
  data <- rbind(
    AGT_row,
    data
  )
}

# 去除可能重复的AGT行
data <- data[
  !duplicated(rownames(data)),
  ,
  drop = FALSE
]

# 将AGT调整到第一行
gene_order <- c(
  target_gene,
  setdiff(rownames(data), target_gene)
)

data <- data[
  gene_order,
  ,
  drop = FALSE
]

# 最终检查
stopifnot(
  target_gene %in% rownames(data),
  rownames(data)[1] == target_gene,
  ncol(data) == ncol(oe_obj)
)

cat("Final input matrix:", nrow(data), "genes x",
    ncol(data), "cells\n")

cat(
  "AGT present before scTenifoldOE:",
  target_gene %in% rownames(data),
  "\n"
)

cat(
  "AGT row position:",
  match(target_gene, rownames(data)),
  "\n"
)

cat(
  "AGT total counts in final matrix:",
  sum(data[target_gene, , drop = FALSE]),
  "\n"
)

############################################################
# 5. AGT虚拟过表达
############################################################

source("refer.scTenifoldOE.R")

set.seed(12345)

result_AGT_OE <- scTenifoldOE(
  countMatrix = data,
  gOE = target_gene,
  oe_factor = 10,
  
  # 关键修改：
  # 前面已经完成人工过滤，因此关闭函数内部QC，
  # 避免scQC再次删除低表达AGT
  qc = FALSE,
  
  nc_nNet = 10,
  nc_nCells = min(
    500,
    ncol(data)
  ),
  nc_nComp = 3,
  nCores = 5
)

saveRDS(
  result_AGT_OE,
  file = "AGT_virtual_OE_Fibroblasts_PD.rds"
)
result_AGT_OE<-readRDS("AGT_virtual_OE_Fibroblasts_PD.rds")
############################################################
# 6. 输出结果
############################################################

df <- result_AGT_OE$diffRegulation

df <- df[
  df$gene != target_gene,
  ,
  drop = FALSE
]

write.csv(
  df,
  file = "AGT_virtual_OE_Fibroblasts_PD_all_results.csv",
  row.names = FALSE
)

sig_df <- df %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  ) %>%
  arrange(p.adj)

write.csv(
  sig_df,
  file = "AGT_virtual_OE_Fibroblasts_PD_significant_genes.csv",
  row.names = FALSE
)

cat("Analysis completed.\n")
cat("Significant genes:", nrow(sig_df), "\n")
############################################################
# 9. 提取差异调控结果
############################################################

df <- result_AGT_OE$diffRegulation

if (is.null(df) || nrow(df) == 0) {
  stop("未获得diffRegulation结果")
}

df <- df[
  df$gene != target_gene,
  ,
  drop = FALSE
]

write.csv(
  df,
  file = "AGT_virtual_OE_Fibroblasts_PD_all_results.csv",
  row.names = FALSE
)

sig_df <- df %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  ) %>%
  arrange(p.adj)

write.csv(
  sig_df,
  file = "AGT_virtual_OE_Fibroblasts_PD_significant_genes.csv",
  row.names = FALSE
)

cat(
  "Significant affected genes:",
  nrow(sig_df),
  "\n"
)

############################################################
# 10. Top 20 affected genes
############################################################

if (nrow(sig_df) > 0) {
  
  top_gene <- sig_df %>%
    arrange(desc(abs(Z))) %>%
    head(20)
  
} else {
  
  warning("没有p.adj < 0.05的基因，Top20图改用全部结果中|Z|最大的20个基因")
  
  top_gene <- df %>%
    filter(!is.na(Z)) %>%
    arrange(desc(abs(Z))) %>%
    head(20)
}

top_gene$gene <- factor(
  top_gene$gene,
  levels = rev(top_gene$gene)
)

p1 <- ggplot(
  top_gene,
  aes(
    x = gene,
    y = Z,
    fill = Z > 0
  )
) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "TRUE" = "#F8766D",
      "FALSE" = "#00BFC4"
    ),
    labels = c(
      "TRUE" = "Positive",
      "FALSE" = "Negative"
    ),
    name = "Z direction"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  labs(
    title = "Genes affected by virtual AGT overexpression",
    subtitle = "PD Fibroblasts",
    x = NULL,
    y = "Perturbation Z-score"
  )

ggsave(
  filename = "AGT_OE_Fibroblasts_PD_top20.pdf",
  plot = p1,
  width = 7,
  height = 5.5
)

ggsave(
  filename = "AGT_OE_Fibroblasts_PD_top20.png",
  plot = p1,
  width = 7,
  height = 5.5,
  dpi = 600
)

############################################################
# 11. 扰动结果散点图
############################################################

df$logFDR <- -log10(
  pmax(df$p.adj, 1e-300)
)

df$status <- ifelse(
  !is.na(df$p.adj) & df$p.adj < 0.05,
  "Significant",
  "NS"
)

label_gene <- df %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.01
  ) %>%
  arrange(p.adj) %>%
  head(15)

p2 <- ggplot(
  df,
  aes(
    x = Z,
    y = logFDR,
    color = status
  )
) +
  geom_point(
    alpha = 0.75,
    size = 1.8
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_text_repel(
    data = label_gene,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Significant" = "#F8766D",
      "NS" = "#00BFC4"
    )
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  labs(
    title = "Virtual AGT overexpression",
    subtitle = "PD Fibroblasts",
    x = "Perturbation Z-score",
    y = "-log10 adjusted P",
    color = NULL
  )

ggsave(
  filename = "AGT_OE_Fibroblasts_PD_volcano.pdf",
  plot = p2,
  width = 7,
  height = 5.5
)

ggsave(
  filename = "AGT_OE_Fibroblasts_PD_volcano.png",
  plot = p2,
  width = 7,
  height = 5.5,
  dpi = 600
)

############################################################
# 12. 输出运行信息
############################################################

cat("\nAGT virtual overexpression analysis completed.\n")
cat("Cell type: Fibroblasts\n")
cat("Group: PD\n")
cat("Target gene: AGT\n")
cat("Overexpression factor: 10\n")
cat("Number of cells:", ncol(data), "\n")
cat("Number of input genes:", nrow(data), "\n")
cat("Significant affected genes:", nrow(sig_df), "\n")