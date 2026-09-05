args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript 03_official_kbet.R <integrated_rds> <output_dir> <local_library>")
}

input_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
local_library <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))

figure_dir <- file.path(output_dir, "results", "figures")
table_dir <- file.path(output_dir, "results", "tables")
log_dir <- file.path(output_dir, "logs")
invisible(lapply(c(figure_dir, table_dir, log_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

log_connection <- file(file.path(log_dir, "official_kbet_run.log"), open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_connection)
}, add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(kBET)
  library(FNN)
  library(ggplot2)
})

set.seed(20260821)
dims_use <- 1:15
n_repeat <- 100
k_requested <- 30L
message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("kBET package version: ", as.character(packageVersion("kBET")))

obj <- readRDS(input_rds)
metadata <- obj[[]]
valid <- which(!is.na(metadata$celltype) & !is.na(metadata$dataset))
strata <- interaction(metadata$celltype[valid], metadata$dataset[valid], drop = TRUE)
stratum_indices <- split(valid, strata)
balanced_indices <- unlist(lapply(stratum_indices, function(index) {
  sample(index, min(length(index), 500L))
}), use.names = FALSE)

run_reduction <- function(embedding, reduction_name) {
  celltypes <- sort(unique(as.character(metadata$celltype[balanced_indices])))
  rows <- lapply(seq_along(celltypes), function(i) {
    celltype_name <- celltypes[[i]]
    index <- balanced_indices[
      as.character(metadata$celltype[balanced_indices]) == celltype_name
    ]
    x <- as.matrix(embedding[index, dims_use, drop = FALSE])
    batch <- factor(as.character(metadata$dataset[index]))
    n <- nrow(x)
    current_k <- min(k_requested, n - 1L)
    expected_counts <- current_k * prop.table(table(batch))

    message(
      reduction_name, " | ", celltype_name,
      " | n=", n, " | k=", current_k
    )
    set.seed(20260821 + i)
    knn <- FNN::get.knn(x, k = current_k, algorithm = "cover_tree")
    result <- kBET::kBET(
      df = x,
      batch = batch,
      k0 = current_k,
      knn = knn,
      do.pca = FALSE,
      heuristic = FALSE,
      n_repeat = n_repeat,
      alpha = 0.05,
      verbose = FALSE,
      plot = FALSE,
      adapt = TRUE
    )
    summary <- result$summary
    data.frame(
      reduction = reduction_name,
      celltype = celltype_name,
      n_balanced = n,
      datasets_present = nlevels(batch),
      k = current_k,
      repeats = n_repeat,
      expected_rejection_rate = summary["mean", "kBET.expected"],
      observed_rejection_rate = summary["mean", "kBET.observed"],
      observed_ci_low = summary["2.5%", "kBET.observed"],
      observed_ci_high = summary["97.5%", "kBET.observed"],
      kbet_signif_fraction = summary["mean", "kBET.signif"],
      composition_limited = any(expected_counts < 5)
    )
  })
  do.call(rbind, rows)
}

official_kbet <- rbind(
  run_reduction(Embeddings(obj, "pca"), "PCA_before"),
  run_reduction(Embeddings(obj, "harmony"), "Harmony_existing"),
  run_reduction(Embeddings(obj, "harmony_explicit"), "Harmony_explicit")
)
write.csv(
  official_kbet,
  file.path(table_dir, "official_kbet_by_celltype.csv"),
  row.names = FALSE
)

official_overall <- do.call(rbind, lapply(split(official_kbet, official_kbet$reduction), function(part) {
  data.frame(
    reduction = part$reduction[[1]],
    weighted_observed_rejection_rate = weighted.mean(
      part$observed_rejection_rate,
      part$n_balanced
    ),
    weighted_expected_rejection_rate = weighted.mean(
      part$expected_rejection_rate,
      part$n_balanced
    ),
    median_celltype_observed_rejection_rate = median(part$observed_rejection_rate),
    celltypes_evaluated = nrow(part),
    composition_limited_celltypes = sum(part$composition_limited)
  )
}))
write.csv(
  official_overall,
  file.path(table_dir, "official_kbet_overall.csv"),
  row.names = FALSE
)

plot_data <- official_kbet
plot_data$reduction <- factor(
  plot_data$reduction,
  levels = c("PCA_before", "Harmony_existing", "Harmony_explicit"),
  labels = c("Before: PCA", "Existing Harmony", "Explicit Harmony")
)
plot_data$celltype <- factor(
  plot_data$celltype,
  levels = rev(unique(plot_data$celltype[order(
    plot_data$observed_rejection_rate[plot_data$reduction == "Explicit Harmony"]
  )]))
)

official_plot <- ggplot(
  plot_data,
  aes(x = observed_rejection_rate, y = celltype, color = reduction)
) +
  geom_errorbarh(
    aes(xmin = observed_ci_low, xmax = observed_ci_high),
    height = 0.15,
    position = position_dodge(width = 0.55)
  ) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c("#B2182B", "#67A9CF", "#2166AC")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Official kBET rejection rate by cell type",
    subtitle = "Lower is better; 100 repeated tests on balanced cell-type strata",
    x = "Observed kBET rejection rate",
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
  file.path(figure_dir, "08_official_kbet_by_celltype.pdf"),
  official_plot,
  width = 11,
  height = 7.5,
  device = cairo_pdf
)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_official_kbet.txt"))
message("Official kBET overall:")
print(official_overall)
message("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
