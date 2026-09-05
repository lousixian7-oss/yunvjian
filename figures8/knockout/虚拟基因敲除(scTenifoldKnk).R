library(Seurat)
library(scTenifoldKnk)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)

set.seed(12345)

#========================
# 1. 读取对象
#========================
yjsl <- readRDS("yjsl_clean_allGSE.rds")

DefaultAssay(yjsl) <- "RNA"

target_gene  <- "AGT"
celltype_use <- "Fibroblasts"
group_use    <- "PD"

# 检查
table(yjsl$cellType, useNA = "ifany")
table(yjsl$group, useNA = "ifany")

#========================
# 2. 提取PD细胞
#========================
ko_obj <- subset(
  yjsl,
  subset =
    cellType == celltype_use &
    group == group_use
)

dim(ko_obj)

#========================
# 3. 提取原始counts
#========================
count_mat <- GetAssayData(
  ko_obj,
  assay = "RNA",
  layer = "counts"
)

if (!target_gene %in% rownames(count_mat)) {
  stop(paste0("目标基因 ", target_gene, " 不在RNA counts矩阵中。"))
}

# 检查目标基因表达
target_detected <- sum(count_mat[target_gene, ] > 0)
target_fraction <- mean(count_mat[target_gene, ] > 0)

cat(
  "阳性细胞数：", target_detected, "\n",
  "阳性比例：", round(target_fraction, 4), "\n"
)


# 非常稀疏，不建议对目标基因本身使用5%阈值
expr_ratio <- Matrix::rowSums(count_mat > 0) / ncol(count_mat)

count_mat_f <- count_mat[
  expr_ratio >= 0.02,
  ,
  drop = FALSE
]

# 无论如何保留
if (!target_gene %in% rownames(count_mat_f)) {
  count_mat_f <- rbind(
    count_mat_f,
    count_mat[target_gene, , drop = FALSE]
  )
}

# 选择HVG
ko_obj <- FindVariableFeatures(
  ko_obj,
  selection.method = "vst",
  nfeatures = 5000,
  verbose = FALSE
)

hvgs <- VariableFeatures(ko_obj)

keep_genes <- unique(
  c(target_gene, hvgs)
)

keep_genes <- intersect(
  keep_genes,
  rownames(count_mat_f)
)

count_mat_f <- count_mat_f[
  keep_genes,
  ,
  drop = FALSE
]

dim(count_mat_f)

nc_cells <- min(
  800,
  ncol(count_mat_f)
)

set.seed(12345)
summary(Matrix::colSums(count_mat))

result_csf2 <- scTenifoldKnk(
  countMatrix = count_mat_f,
  gKO = target_gene,
  qc = TRUE,
  qc_mtThreshold = 0.1,
  qc_minLSize = 1000,
  nc_nNet = 10,
  nc_nCells = nc_cells
)

saveRDS(
  result_csf2,
  "AGT_virtual_KO_Endothelial_PD.rds"
)

df_ko <- as.data.frame(
  result_csf2$diffRegulation
)

df_ko <- df_ko[
  df_ko$gene != target_gene,
  ,
  drop = FALSE
]

colnames(df_ko)
head(df_ko)

sig_ko <- df_ko %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  ) %>%
  arrange(
    p.adj,
    desc(abs(Z))
  )

nrow(sig_ko)

write.csv(
  df_ko,
  "CSF2_virtual_KO_all_results.csv",
  row.names = FALSE
)

write.csv(
  sig_ko,
  "CSF2_virtual_KO_significant_genes.csv",
  row.names = FALSE
)

top_n <- min(
  15,
  nrow(sig_ko)
)

if (top_n > 0) {
  
  top_genes <- sig_ko %>%
    arrange(
      desc(distance)
    ) %>%
    slice_head(
      n = top_n
    )
  
  top_genes$gene <- factor(
    top_genes$gene,
    levels = rev(top_genes$gene)
  )
  
  p_bar <- ggplot(
    top_genes,
    aes(
      x = gene,
      y = distance
    )
  ) +
    geom_col(
      width = 0.72
    ) +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      x = NULL,
      y = "Manifold distance",
      title = "Top genes affected by virtual CSF2 knockout"
    )
  
  p_bar
  
  ggsave(
    "CSF2_virtual_KO_top_genes.pdf",
    p_bar,
    width = 6,
    height = 5
  )
}

df_ko <- df_ko %>%
  mutate(
    log10_FDR = -log10(
      pmax(p.adj, 1e-300)
    ),
    Significance = ifelse(
      p.adj < 0.05 & abs(Z) >= 2,
      "Significant",
      "Not significant"
    )
  )

label_genes <- df_ko %>%
  filter(
    p.adj < 0.01,
    abs(Z) >= 2
  ) %>%
  arrange(
    p.adj
  ) %>%
  slice_head(
    n = 15
  )

p_perturb <- ggplot(
  df_ko,
  aes(
    x = Z,
    y = log10_FDR,
    color = Significance
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.2
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = c(-2, 2),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Significant" = "#F8766D",
      "Not significant" = "#00BFC4"
    )
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Perturbation Z-score",
    y = expression(-log[10]("adjusted P value")),
    title = "Predicted transcriptional effects of virtual CSF2 knockout"
  ) +
  theme(
    legend.position = "none"
  )

p_perturb

ggsave(
  "AGT_virtual_KO_perturbation_plot.pdf",
  p_perturb,
  width = 6,
  height = 5
)