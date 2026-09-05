args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    "Usage: Rscript 01_batch_integration_revision.R ",
    "<input_rds> <output_dir> <local_library>"
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

log_file <- file.path(log_dir, "batch_integration_run.log")
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
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)
dims_use <- 1:15
dataset_colors <- c(
  GSE152042 = "#F8766D",
  GSE164241 = "#00BA38",
  GSE171213 = "#619CFF"
)

message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("Input: ", input_rds)
message("Output: ", output_dir)

obj <- readRDS(input_rds)
required_metadata <- c("dataset", "sample", "group", "celltype")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata) > 0L) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}
missing_reductions <- setdiff(c("pca", "harmony", "umap"), Reductions(obj))
if (length(missing_reductions) > 0L) {
  stop("Missing reductions: ", paste(missing_reductions, collapse = ", "))
}
if (max(dims_use) > ncol(Embeddings(obj, "pca")) ||
    max(dims_use) > ncol(Embeddings(obj, "harmony"))) {
  stop("Requested dimensions are not available in PCA or Harmony.")
}

message("Object: ", nrow(obj), " genes x ", ncol(obj), " cells")
message("Reductions: ", paste(Reductions(obj), collapse = ", "))

# Before integration: PCA-based UMAP.
obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = dims_use,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = 42,
  reduction.name = "umap_before",
  reduction.key = "UMAPbefore_",
  verbose = TRUE
)

# Existing Harmony embedding, regenerated with named UMAP parameters.
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = dims_use,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = 42,
  reduction.name = "umap_harmony_existing",
  reduction.key = "UMAPexisting_",
  verbose = TRUE
)

# Explicit Harmony rerun. This preserves the original `harmony` reduction.
convergence_pdf <- file.path(figure_dir, "01_harmony_explicit_convergence.pdf")
grDevices::pdf(convergence_pdf, width = 7, height = 5, useDingbats = FALSE)
obj <- tryCatch(
  RunHarmony(
    object = obj,
    group.by.vars = "dataset",
    reduction.use = "pca",
    dims.use = dims_use,
    theta = 2,
    lambda = 1,
    sigma = 0.1,
    nclust = NULL,
    max_iter = 10,
    early_stop = TRUE,
    ncores = 1,
    plot_convergence = TRUE,
    verbose = TRUE,
    reduction.save = "harmony_explicit",
    project.dim = TRUE,
    .options = harmony_options(
      alpha = 0.2,
      tau = 0,
      block.size = 0.05,
      max.iter.cluster = 4,
      epsilon.cluster = 1e-3,
      epsilon.harmony = 1e-2,
      batch.prop.cutoff = 1e-5
    )
  ),
  error = function(error) {
    grDevices::dev.off()
    stop(error)
  }
)
grDevices::dev.off()

obj <- RunUMAP(
  obj,
  reduction = "harmony_explicit",
  dims = dims_use,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = 42,
  reduction.name = "umap_harmony_explicit",
  reduction.key = "UMAPexplicit_",
  verbose = TRUE
)

plot_dataset <- function(reduction, title) {
  DimPlot(
    obj,
    reduction = reduction,
    group.by = "dataset",
    cols = dataset_colors,
    raster = TRUE,
    shuffle = TRUE,
    seed = 42
  ) +
    ggtitle(title) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

plot_sample <- function(reduction, title) {
  DimPlot(
    obj,
    reduction = reduction,
    group.by = "sample",
    raster = TRUE,
    shuffle = TRUE,
    seed = 42
  ) +
    ggtitle(title) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend()
}

plot_celltype <- function(reduction, title) {
  DimPlot(
    obj,
    reduction = reduction,
    group.by = "celltype",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle(title) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

dataset_comparison <-
  plot_dataset("umap_before", "Before integration: PCA") |
  plot_dataset("umap_harmony_existing", "After integration: existing Harmony") |
  plot_dataset("umap_harmony_explicit", "After integration: explicit Harmony rerun")

sample_comparison <-
  plot_sample("umap_before", "Before integration: sample") |
  plot_sample("umap_harmony_explicit", "After integration: sample")

celltype_comparison <-
  plot_celltype("umap_before", "Before integration: cell type") |
  plot_celltype("umap_harmony_explicit", "After integration: cell type")

split_before <- DimPlot(
  obj,
  reduction = "umap_before",
  group.by = "dataset",
  split.by = "dataset",
  cols = dataset_colors,
  raster = TRUE,
  ncol = 3
) + plot_annotation(title = "Before integration: dataset-split UMAP")

split_after <- DimPlot(
  obj,
  reduction = "umap_harmony_explicit",
  group.by = "dataset",
  split.by = "dataset",
  cols = dataset_colors,
  raster = TRUE,
  ncol = 3
) + plot_annotation(title = "After integration: dataset-split UMAP")

ggsave(
  file.path(figure_dir, "02_dataset_before_existing_explicit.pdf"),
  dataset_comparison,
  width = 18,
  height = 6,
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "03_sample_before_after.pdf"),
  sample_comparison,
  width = 14,
  height = 6,
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "04_celltype_before_after.pdf"),
  celltype_comparison,
  width = 16,
  height = 7,
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "05_dataset_split_before_after.pdf"),
  split_before / split_after,
  width = 16,
  height = 10,
  device = cairo_pdf
)

# Composition tables.
metadata <- obj[[]]
dataset_group_counts <- as.data.frame(with(metadata, table(dataset, group)))
names(dataset_group_counts) <- c("dataset", "group", "cell_count")
write.csv(
  dataset_group_counts,
  file.path(table_dir, "dataset_group_cell_counts.csv"),
  row.names = FALSE
)

dataset_sample_group_counts <- as.data.frame(
  with(metadata, table(dataset, sample, group))
)
names(dataset_sample_group_counts) <- c("dataset", "sample", "group", "cell_count")
dataset_sample_group_counts <- subset(dataset_sample_group_counts, cell_count > 0)
write.csv(
  dataset_sample_group_counts,
  file.path(table_dir, "dataset_sample_group_cell_counts.csv"),
  row.names = FALSE
)

celltype_for_table <- ifelse(
  is.na(metadata$celltype),
  "Unannotated",
  as.character(metadata$celltype)
)
celltype_dataset_counts <- as.data.frame(
  table(celltype = celltype_for_table, dataset = metadata$dataset)
)
names(celltype_dataset_counts) <- c("celltype", "dataset", "cell_count")
celltype_dataset_counts$celltype_total <- ave(
  celltype_dataset_counts$cell_count,
  celltype_dataset_counts$celltype,
  FUN = sum
)
celltype_dataset_counts$percent_within_celltype <-
  100 * celltype_dataset_counts$cell_count / celltype_dataset_counts$celltype_total
write.csv(
  celltype_dataset_counts,
  file.path(table_dir, "celltype_dataset_composition.csv"),
  row.names = FALSE
)

summarize_qc <- function(values) {
  c(
    median = median(values, na.rm = TRUE),
    q05 = unname(quantile(values, 0.05, na.rm = TRUE)),
    q95 = unname(quantile(values, 0.95, na.rm = TRUE)),
    max = max(values, na.rm = TRUE)
  )
}

qc_rows <- lapply(split(metadata, metadata$dataset), function(part) {
  data.frame(
    dataset = part$dataset[[1]],
    cells_postfilter = nrow(part),
    percent_mt_median = median(part$percent.mt, na.rm = TRUE),
    percent_mt_q95 = unname(quantile(part$percent.mt, 0.95, na.rm = TRUE)),
    percent_mt_max = max(part$percent.mt, na.rm = TRUE),
    nFeature_RNA_median = median(part$nFeature_RNA, na.rm = TRUE),
    nFeature_RNA_q95 = unname(quantile(part$nFeature_RNA, 0.95, na.rm = TRUE)),
    nCount_RNA_median = median(part$nCount_RNA, na.rm = TRUE),
    nCount_RNA_q95 = unname(quantile(part$nCount_RNA, 0.95, na.rm = TRUE))
  )
})
qc_summary <- do.call(rbind, qc_rows)
write.csv(
  qc_summary,
  file.path(table_dir, "postfilter_qc_summary_by_dataset.csv"),
  row.names = FALSE
)

# Save the completed reductions before downstream metric calculations so that
# a reporting-only failure cannot discard the expensive Harmony/UMAP results.
saveRDS(
  obj,
  file.path(object_dir, "yjsl_batch_integration_revision.rds"),
  compress = "gzip"
)

# Balanced, cell-type-stratified batch metrics.
set.seed(20260821)
valid <- which(!is.na(metadata$celltype) & !is.na(metadata$dataset))
strata <- interaction(
  metadata$celltype[valid],
  metadata$dataset[valid],
  drop = TRUE
)
stratum_indices <- split(valid, strata)
balanced_indices <- unlist(lapply(stratum_indices, function(index) {
  sample(index, min(length(index), 500L))
}), use.names = FALSE)

calculate_metrics <- function(embedding, reduction_name, indices, metadata, k = 30L) {
  celltypes <- sort(unique(as.character(metadata$celltype[indices])))
  rows <- lapply(celltypes, function(celltype_name) {
    index <- indices[metadata$celltype[indices] == celltype_name]
    x <- as.matrix(embedding[index, dims_use, drop = FALSE])
    labels <- as.character(metadata$dataset[index])
    n <- nrow(x)
    current_k <- min(k, n - 1L)

    distances <- as.matrix(stats::dist(x, method = "euclidean"))
    diag(distances) <- Inf
    neighbors <- vapply(seq_len(n), function(i) {
      order(distances[i, ], method = "radix")[seq_len(current_k)]
    }, integer(current_k))

    neighbor_labels <- matrix(labels[neighbors], nrow = current_k, ncol = n)
    cross30 <- mean(neighbor_labels != rep(labels, each = current_k))

    dataset_levels <- sort(unique(labels))
    ilisi_per_cell <- vapply(seq_len(n), function(i) {
      proportions <- table(factor(neighbor_labels[, i], levels = dataset_levels)) /
        current_k
      1 / sum(proportions^2)
    }, numeric(1))

    silhouette_per_cell <- vapply(seq_len(n), function(i) {
      same <- which(labels == labels[[i]])
      same <- setdiff(same, i)
      a <- if (length(same) > 0L) mean(distances[i, same]) else 0
      other_means <- vapply(
        setdiff(dataset_levels, labels[[i]]),
        function(other) mean(distances[i, labels == other]),
        numeric(1)
      )
      b <- min(other_means)
      if (max(a, b) == 0) 0 else (b - a) / max(a, b)
    }, numeric(1))

    proportions <- prop.table(table(labels))
    expected_cross <- 1 - sum(proportions^2)

    data.frame(
      reduction = reduction_name,
      celltype = celltype_name,
      n_balanced = n,
      expected_cross = expected_cross,
      cross_dataset_30nn = cross30,
      cross_to_expected_ratio = cross30 / expected_cross,
      iLISI = mean(ilisi_per_cell),
      dataset_silhouette = mean(silhouette_per_cell)
    )
  })
  do.call(rbind, rows)
}

metric_results <- rbind(
  calculate_metrics(Embeddings(obj, "pca"), "PCA_before", balanced_indices, metadata),
  calculate_metrics(Embeddings(obj, "harmony"), "Harmony_existing", balanced_indices, metadata),
  calculate_metrics(
    Embeddings(obj, "harmony_explicit"),
    "Harmony_explicit",
    balanced_indices,
    metadata
  )
)
write.csv(
  metric_results,
  file.path(table_dir, "batch_metrics_by_celltype.csv"),
  row.names = FALSE
)

numeric_metrics <- c(
  "expected_cross",
  "cross_dataset_30nn",
  "cross_to_expected_ratio",
  "iLISI",
  "dataset_silhouette"
)
overall_metrics <- do.call(rbind, lapply(split(metric_results, metric_results$reduction), function(part) {
  output <- data.frame(reduction = part$reduction[[1]])
  for (metric in numeric_metrics) {
    output[[metric]] <- weighted.mean(part[[metric]], part$n_balanced)
  }
  output
}))
write.csv(
  overall_metrics,
  file.path(table_dir, "batch_metrics_overall.csv"),
  row.names = FALSE
)

existing_harmony <- Embeddings(obj, "harmony")[, dims_use, drop = FALSE]
explicit_harmony <- Embeddings(obj, "harmony_explicit")[, dims_use, drop = FALSE]
harmony_comparison <- data.frame(
  dimension = dims_use,
  pearson_correlation = vapply(dims_use, function(dimension) {
    cor(existing_harmony[, dimension], explicit_harmony[, dimension])
  }, numeric(1)),
  root_mean_squared_difference = vapply(dims_use, function(dimension) {
    sqrt(mean((existing_harmony[, dimension] - explicit_harmony[, dimension])^2))
  }, numeric(1))
)
write.csv(
  harmony_comparison,
  file.path(table_dir, "existing_vs_explicit_harmony_dimensions.csv"),
  row.names = FALSE
)

parameters <- data.frame(
  parameter = c(
    "seed", "group.by.vars", "reduction", "dims.use", "theta", "lambda",
    "sigma", "max_iter", "early_stop", "ncores", "alpha", "tau",
    "block.size", "max.iter.cluster", "epsilon.cluster", "epsilon.harmony",
    "batch.prop.cutoff", "UMAP n.neighbors", "UMAP min.dist", "UMAP metric"
  ),
  value = c(
    "42", "dataset", "pca", "1:15", "2", "1", "0.1", "10", "TRUE",
    "1", "0.2", "0", "0.05", "4", "1e-3", "1e-2", "1e-5", "30",
    "0.3", "cosine"
  )
)
write.csv(
  parameters,
  file.path(table_dir, "integration_parameters.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(log_dir, "sessionInfo.txt")
)

message("Overall metrics:")
print(overall_metrics)
message("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
