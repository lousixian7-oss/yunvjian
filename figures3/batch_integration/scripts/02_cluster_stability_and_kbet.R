args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    "Usage: Rscript 02_cluster_stability_and_kbet.R ",
    "<integrated_rds> <output_dir> <local_library>"
  )
}

input_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
local_library <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))

figure_dir <- file.path(output_dir, "results", "figures")
table_dir <- file.path(output_dir, "results", "tables")
object_dir <- file.path(output_dir, "results", "objects")
log_dir <- file.path(output_dir, "logs")
invisible(lapply(
  c(figure_dir, table_dir, object_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

log_file <- file.path(log_dir, "cluster_stability_kbet_run.log")
log_connection <- file(log_file, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_connection)
}, add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)
dims_use <- 1:15
message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("Input: ", input_rds)

obj <- readRDS(input_rds)
metadata <- obj[[]]
required_metadata <- c("dataset", "sample", "celltype", "seurat_clusters")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata) > 0L) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}
if (!"harmony_explicit" %in% Reductions(obj)) {
  stop("Reduction 'harmony_explicit' is missing.")
}

obj$clusters_original <- factor(as.character(obj$seurat_clusters))
obj <- FindNeighbors(
  obj,
  reduction = "harmony_explicit",
  dims = dims_use,
  k.param = 20,
  graph.name = c("explicit_harmony_nn", "explicit_harmony_snn"),
  verbose = TRUE
)
obj <- FindClusters(
  obj,
  graph.name = "explicit_harmony_snn",
  resolution = 0.3,
  algorithm = 1,
  random.seed = 42,
  cluster.name = "clusters_explicit",
  verbose = TRUE
)

old_cluster <- factor(as.character(obj$clusters_original))
new_cluster <- factor(as.character(obj$clusters_explicit))

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(sum(tab))
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- (sum_rows + sum_cols) / 2
  if (maximum == expected) return(1)
  (sum_cells - expected) / (maximum - expected)
}

normalized_mutual_information <- function(x, y) {
  tab <- table(x, y)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(vapply(seq_len(nrow(nz)), function(i) {
    row <- nz[i, 1]
    col <- nz[i, 2]
    pxy[row, col] * log(pxy[row, col] / (px[row] * py[col]))
  }, numeric(1)))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(1)
  mi / sqrt(hx * hy)
}

overlap <- table(original_cluster = old_cluster, explicit_cluster = new_cluster)
old_to_new_purity <- sum(apply(overlap, 1, max)) / sum(overlap)
new_to_old_purity <- sum(apply(overlap, 2, max)) / sum(overlap)

celltype <- factor(ifelse(is.na(obj$celltype), "Unannotated", as.character(obj$celltype)))
celltype_by_new <- table(celltype = celltype, explicit_cluster = new_cluster)
new_cluster_celltype_purity <- sum(apply(celltype_by_new, 2, max)) /
  sum(celltype_by_new)

new_to_old_map <- apply(overlap, 2, function(values) {
  rownames(overlap)[which.max(values)]
})
mapped_new <- factor(
  unname(new_to_old_map[as.character(new_cluster)]),
  levels = levels(old_cluster)
)
cell_mapping_accuracy <- mean(mapped_new == old_cluster)

sample_levels <- sort(unique(as.character(obj$sample)))
cluster_levels <- levels(old_cluster)
sample_profile_rows <- lapply(sample_levels, function(sample_name) {
  index <- which(as.character(obj$sample) == sample_name)
  old_profile <- prop.table(table(factor(old_cluster[index], levels = cluster_levels)))
  new_profile <- prop.table(table(factor(mapped_new[index], levels = cluster_levels)))
  correlation <- suppressWarnings(cor(as.numeric(old_profile), as.numeric(new_profile)))
  if (!is.finite(correlation)) correlation <- NA_real_
  data.frame(
    sample = sample_name,
    n_cells = length(index),
    pearson_cluster_profile = correlation,
    total_variation_distance = 0.5 * sum(abs(old_profile - new_profile))
  )
})
sample_profile <- do.call(rbind, sample_profile_rows)

stability_summary <- data.frame(
  metric = c(
    "cells", "original_cluster_count", "explicit_cluster_count",
    "adjusted_rand_index", "normalized_mutual_information",
    "original_to_explicit_overlap_purity",
    "explicit_to_original_overlap_purity",
    "cell_mapping_accuracy_after_majority_map",
    "new_cluster_celltype_purity",
    "median_sample_cluster_profile_correlation",
    "median_sample_total_variation_distance"
  ),
  value = c(
    ncol(obj), nlevels(old_cluster), nlevels(new_cluster),
    adjusted_rand_index(old_cluster, new_cluster),
    normalized_mutual_information(old_cluster, new_cluster),
    old_to_new_purity, new_to_old_purity, cell_mapping_accuracy,
    new_cluster_celltype_purity,
    median(sample_profile$pearson_cluster_profile, na.rm = TRUE),
    median(sample_profile$total_variation_distance, na.rm = TRUE)
  )
)

write.csv(stability_summary, file.path(table_dir, "cluster_stability_summary.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(overlap), file.path(table_dir, "cluster_overlap_counts.csv"))
write.csv(
  as.data.frame.matrix(prop.table(overlap, margin = 1)),
  file.path(table_dir, "cluster_overlap_row_proportions.csv")
)
write.csv(
  as.data.frame.matrix(celltype_by_new),
  file.path(table_dir, "new_cluster_celltype_counts.csv")
)
write.csv(sample_profile, file.path(table_dir, "sample_cluster_profile_stability.csv"), row.names = FALSE)

plot_old <- DimPlot(
  obj,
  reduction = "umap_harmony_explicit",
  group.by = "clusters_original",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("Original clustering on explicit Harmony UMAP") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

plot_new <- DimPlot(
  obj,
  reduction = "umap_harmony_explicit",
  group.by = "clusters_explicit",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("Recomputed graph/clusters on explicit Harmony") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(
  file.path(figure_dir, "06_original_vs_explicit_clusters.pdf"),
  plot_old | plot_new,
  width = 15,
  height = 7,
  device = cairo_pdf
)

overlap_long <- as.data.frame(prop.table(overlap, margin = 1))
names(overlap_long) <- c("original_cluster", "explicit_cluster", "row_fraction")
overlap_heatmap <- ggplot(
  overlap_long,
  aes(x = explicit_cluster, y = original_cluster, fill = row_fraction)
) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1)) +
  labs(
    title = "Cluster overlap after explicit Harmony graph reconstruction",
    subtitle = "Each row sums to 1; darker cells indicate stronger one-to-one preservation",
    x = "Explicit Harmony cluster",
    y = "Original cluster",
    fill = "Row fraction"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(figure_dir, "07_cluster_overlap_heatmap.pdf"),
  overlap_heatmap,
  width = 10,
  height = 8,
  device = cairo_pdf
)

# kBET-style test: within each annotated cell type, compare the local 30-NN
# dataset composition against that cell type's global dataset composition.
# Pearson chi-square statistics are calibrated by a multinomial Monte Carlo
# null, which remains valid when an expected local batch count is below five.
set.seed(20260821)
valid <- which(!is.na(obj$celltype) & !is.na(obj$dataset))
strata <- interaction(obj$celltype[valid], obj$dataset[valid], drop = TRUE)
stratum_indices <- split(valid, strata)
balanced_indices <- unlist(lapply(stratum_indices, function(index) {
  sample(index, min(length(index), 500L))
}), use.names = FALSE)

kbet_style_one <- function(embedding, reduction_name, indices, metadata, k = 30L) {
  celltypes <- sort(unique(as.character(metadata$celltype[indices])))
  rows <- lapply(seq_along(celltypes), function(celltype_number) {
    celltype_name <- celltypes[[celltype_number]]
    index <- indices[as.character(metadata$celltype[indices]) == celltype_name]
    x <- as.matrix(embedding[index, dims_use, drop = FALSE])
    labels <- as.character(metadata$dataset[index])
    n <- nrow(x)
    current_k <- min(k, n - 1L)
    levels_batch <- sort(unique(labels))
    batch_prob <- prop.table(table(factor(labels, levels = levels_batch)))

    distances <- as.matrix(stats::dist(x, method = "euclidean"))
    diag(distances) <- Inf
    neighbors <- vapply(seq_len(n), function(i) {
      order(distances[i, ], method = "radix")[seq_len(current_k)]
    }, integer(current_k))
    neighbor_labels <- matrix(labels[neighbors], nrow = current_k, ncol = n)
    observed <- vapply(seq_len(n), function(i) {
      as.numeric(table(factor(neighbor_labels[, i], levels = levels_batch)))
    }, numeric(length(levels_batch)))
    expected <- current_k * as.numeric(batch_prob)
    observed_stat <- colSums((observed - expected)^2 / expected)

    set.seed(20260821 + celltype_number)
    null_counts <- t(rmultinom(100000L, size = current_k, prob = batch_prob))
    null_stat <- rowSums((null_counts - matrix(
      expected,
      nrow = nrow(null_counts),
      ncol = length(expected),
      byrow = TRUE
    ))^2 / matrix(
      expected,
      nrow = nrow(null_counts),
      ncol = length(expected),
      byrow = TRUE
    ))
    p_values <- (vapply(observed_stat, function(value) {
      sum(null_stat >= value)
    }, integer(1)) + 1) / (length(null_stat) + 1)
    rejected <- p_values < 0.05

    set.seed(20260821 + 1000L + celltype_number)
    test_size <- min(n, max(25L, ceiling(0.1 * n)))
    repeat_rates <- replicate(100L, mean(sample(rejected, test_size, replace = FALSE)))

    data.frame(
      reduction = reduction_name,
      celltype = celltype_name,
      n_balanced = n,
      datasets_present = length(levels_batch),
      k = current_k,
      tested_per_repeat = test_size,
      rejection_rate = mean(repeat_rates),
      ci_low = unname(quantile(repeat_rates, 0.025)),
      ci_high = unname(quantile(repeat_rates, 0.975)),
      median_local_p = median(p_values),
      composition_limited = any(expected < 5)
    )
  })
  do.call(rbind, rows)
}

metadata <- obj[[]]
kbet_style <- rbind(
  kbet_style_one(Embeddings(obj, "pca"), "PCA_before", balanced_indices, metadata),
  kbet_style_one(Embeddings(obj, "harmony"), "Harmony_existing", balanced_indices, metadata),
  kbet_style_one(
    Embeddings(obj, "harmony_explicit"),
    "Harmony_explicit",
    balanced_indices,
    metadata
  )
)
write.csv(kbet_style, file.path(table_dir, "kbet_style_by_celltype.csv"), row.names = FALSE)

kbet_overall <- do.call(rbind, lapply(split(kbet_style, kbet_style$reduction), function(part) {
  data.frame(
    reduction = part$reduction[[1]],
    weighted_rejection_rate = weighted.mean(part$rejection_rate, part$n_balanced),
    median_celltype_rejection_rate = median(part$rejection_rate),
    celltypes_evaluated = nrow(part),
    composition_limited_celltypes = sum(part$composition_limited)
  )
}))
write.csv(kbet_overall, file.path(table_dir, "kbet_style_overall.csv"), row.names = FALSE)

kbet_plot_data <- kbet_style
kbet_plot_data$reduction <- factor(
  kbet_plot_data$reduction,
  levels = c("PCA_before", "Harmony_existing", "Harmony_explicit"),
  labels = c("Before: PCA", "Existing Harmony", "Explicit Harmony")
)
kbet_plot <- ggplot(
  kbet_plot_data,
  aes(x = rejection_rate, y = reorder(celltype, rejection_rate), color = reduction)
) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey45") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c("#B2182B", "#67A9CF", "#2166AC")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Cell-type-stratified kBET-style rejection rate",
    subtitle = "Lower is better; dashed line marks the nominal 0.05 rejection level",
    x = "Local batch-composition rejection rate",
    y = "Cell type",
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  file.path(figure_dir, "08_kbet_style_by_celltype.pdf"),
  kbet_plot,
  width = 11,
  height = 7.5,
  device = cairo_pdf
)

decision <- data.frame(
  rule = c(
    "ARI >= 0.90",
    "cell mapping accuracy >= 0.95",
    "new cluster cell-type purity >= 0.95",
    "median sample cluster-profile correlation >= 0.95"
  ),
  observed = c(
    stability_summary$value[stability_summary$metric == "adjusted_rand_index"],
    cell_mapping_accuracy,
    new_cluster_celltype_purity,
    median(sample_profile$pearson_cluster_profile, na.rm = TRUE)
  ),
  passes = c(
    stability_summary$value[stability_summary$metric == "adjusted_rand_index"] >= 0.90,
    cell_mapping_accuracy >= 0.95,
    new_cluster_celltype_purity >= 0.95,
    median(sample_profile$pearson_cluster_profile, na.rm = TRUE) >= 0.95
  )
)
write.csv(decision, file.path(table_dir, "stability_decision_rules.csv"), row.names = FALSE)

saveRDS(
  obj,
  file.path(object_dir, "yjsl_batch_integration_revision_clustered.rds"),
  compress = "gzip"
)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_stability.txt"))
message("Stability summary:")
print(stability_summary)
message("kBET-style overall:")
print(kbet_overall)
message("Decision rules:")
print(decision)
message("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
