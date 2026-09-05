#!/usr/bin/env Rscript

options(width = 220)
.libPaths(c("G:/gurobi/Rlib45", .libPaths()))

suppressPackageStartupMessages({
  library(mistyR)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

set.seed(20260824)
theme_set(theme_classic(base_size = 13))

spatial_file <- "stRNA_RCTD_refinedFib_ECM_PD.rds"
weights_file <- "RCTD_PD_refinedFib_ECM_weights_full.csv"
results_dir <- "MISTy_PD_RCTD_refinedFib_ECM_full"
log_file <- "MISTy_PD_RCTD_refinedFib_ECM_full.log"

stopifnot(file.exists(spatial_file), file.exists(weights_file))

sink(log_file, split = TRUE)
on.exit({ while (sink.number() > 0L) sink() }, add = TRUE)

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("mistyR:", as.character(packageVersion("mistyR")), "\n")
cat("Input weights: weights_full (continuous RCTD weights)\n")

stRNA <- readRDS(spatial_file)
cat("Spatial object class:", paste(class(stRNA), collapse = ", "), "\n")

weights_df <- read.csv(weights_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"spot_id" %in% colnames(weights_df)) stop("weights_full CSV lacks spot_id")
rownames(weights_df) <- weights_df$spot_id
weights_df$spot_id <- NULL
rctd_mat <- as.matrix(weights_df)
storage.mode(rctd_mat) <- "double"

if (any(!is.finite(rctd_mat)) || any(rctd_mat < 0)) stop("Invalid RCTD weights")
keep_spots <- rowSums(rctd_mat) > 0
rctd_mat <- rctd_mat[keep_spots, , drop = FALSE]
rctd_mat <- rctd_mat / rowSums(rctd_mat)
keep_types <- apply(rctd_mat, 2, sd) > 0
rctd_mat <- rctd_mat[, keep_types, drop = FALSE]

if (!all(c("ECM_Fib", "Activated_Fib") %in% colnames(rctd_mat))) {
  stop("ECM_Fib or Activated_Fib is absent from weights_full")
}

# Extract coordinates without requiring Seurat to be installed.
images <- attr(stRNA, "images", exact = TRUE)
if (!length(images)) stop("Spatial object has no image")
image_attrs <- attributes(images[[1]])
centroid_attrs <- attributes(image_attrs$boundaries[["centroids"]])
coords <- as.data.frame(centroid_attrs$coords)
rownames(coords) <- centroid_attrs$cells
colnames(coords)[1:2] <- c("x", "y")
coords <- coords[, c("x", "y"), drop = FALSE]

common_spots <- intersect(rownames(rctd_mat), rownames(coords))
if (!length(common_spots)) stop("No common spots between RCTD and spatial coordinates")
rctd_mat <- rctd_mat[common_spots, , drop = FALSE]
coords <- coords[common_spots, , drop = FALSE]
stopifnot(identical(rownames(rctd_mat), rownames(coords)))

name_map <- data.frame(
  old_name = colnames(rctd_mat),
  new_name = make.names(colnames(rctd_mat), unique = TRUE),
  stringsAsFactors = FALSE
)
colnames(rctd_mat) <- name_map$new_name
write.csv(name_map, "MISTy_PD_celltype_name_map_refinedFib_ECM.csv", row.names = FALSE)

input_summary <- data.frame(
  celltype = colnames(rctd_mat),
  mean_proportion = colMeans(rctd_mat),
  sd_proportion = apply(rctd_mat, 2, sd),
  min_proportion = apply(rctd_mat, 2, min),
  max_proportion = apply(rctd_mat, 2, max),
  row.names = NULL
)
write.csv(input_summary, "MISTy_PD_RCTD_refinedFib_ECM_input_summary.csv", row.names = FALSE)

cat("Aligned spots:", nrow(rctd_mat), "\n")
cat("Cell types:", ncol(rctd_mat), "\n")
cat("Nearest-neighbor threshold: 50 coordinate units\n")
cat("Paraview length scale: 100 coordinate units\n")

misty_views <- create_initial_view(as.data.frame(rctd_mat)) %>%
  add_juxtaview(
    positions = coords,
    neighbor.thr = 50,
    prefix = "juxta"
  ) %>%
  add_paraview(
    positions = coords,
    l = 100,
    prefix = "para"
  )

if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::multisession, workers = 4)
  on.exit(future::plan(future::sequential), add = TRUE)
}

misty_result <- run_misty(
  misty_views,
  results.folder = results_dir,
  seed = 20260824,
  cv.folds = 10,
  cached = FALSE
)
saveRDS(misty_result, "MISTy_PD_RCTD_refinedFib_ECM_run_return.rds")

misty_collect <- collect_results(results_dir)
saveRDS(misty_collect, "MISTy_PD_RCTD_refinedFib_ECM_collect.rds")
write.csv(misty_collect$improvements, "MISTy_PD_improvements.csv", row.names = FALSE)
write.csv(misty_collect$contributions, "MISTy_PD_contributions.csv", row.names = FALSE)
write.csv(
  misty_collect$importances.aggregated,
  "MISTy_PD_importances_aggregated.csv",
  row.names = FALSE
)

save_misty_plot <- function(filename, expression, width, height) {
  grDevices::pdf(filename, width = width, height = height, onefile = TRUE)
  tryCatch(
    force(expression),
    finally = grDevices::dev.off()
  )
  invisible(NULL)
}

# mistyR 1.18 draws these plots as a side effect and returns a tibble, so use
# an explicit PDF device instead of passing the return value to ggsave().
save_misty_plot(
  "MISTy_PD_RCTD_refinedFib_ECM_improvement.pdf",
  plot_improvement_stats(misty_collect),
  width = 7.2,
  height = 5.2
)

view_names <- unique(misty_collect$importances.aggregated$view)
for (view_name in view_names) {
  safe_view <- gsub("[^A-Za-z0-9_.-]", "_", view_name)
  tryCatch(
    save_misty_plot(
      paste0("MISTy_PD_RCTD_refinedFib_ECM_interactions_", safe_view, ".pdf"),
      plot_interaction_heatmap(misty_collect, view = view_name, cutoff = 0),
      width = 7.2,
      height = 6.4
    ),
    error = function(e) {
      message("Skipping interaction heatmap for ", view_name, ": ", conditionMessage(e))
    }
  )
}

plot_df <- cbind(coords, as.data.frame(rctd_mat))
focus_long <- plot_df %>%
  mutate(spot_id = rownames(plot_df)) %>%
  select(spot_id, x, y, ECM_Fib, Activated_Fib) %>%
  pivot_longer(c(ECM_Fib, Activated_Fib), names_to = "celltype", values_to = "proportion")

p_focus <- ggplot(focus_long, aes(x = x, y = y, color = proportion)) +
  geom_point(size = 1.2, alpha = 0.95) +
  scale_color_viridis_c(option = "turbo") +
  scale_y_reverse() +
  coord_fixed() +
  facet_wrap(~celltype, nrow = 1) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", strip.text = element_text(face = "bold")) +
  labs(color = "RCTD proportion")
ggsave("PD_ECM_Fib_Activated_Fib_spatial_MISTy_input.pdf", p_focus, width = 8.8, height = 4.7)

cor_df <- as.data.frame(rctd_mat)[, c("ECM_Fib", "Activated_Fib"), drop = FALSE]
pearson <- cor.test(cor_df$ECM_Fib, cor_df$Activated_Fib, method = "pearson")
spearman <- suppressWarnings(
  cor.test(cor_df$ECM_Fib, cor_df$Activated_Fib, method = "spearman", exact = FALSE)
)
cor_stats <- data.frame(
  n_spots = nrow(cor_df),
  pearson_r = unname(pearson$estimate),
  pearson_p = pearson$p.value,
  spearman_rho = unname(spearman$estimate),
  spearman_p = spearman$p.value
)
write.csv(cor_stats, "PD_ECM_Fib_Activated_Fib_spot_correlation_stats.csv", row.names = FALSE)

format_p <- function(p) {
  if (p < 2.2e-16) "< 2.2e-16" else format.pval(p, digits = 3)
}
cor_label <- sprintf(
  "n = %d\nPearson r = %.3f, p %s\nSpearman rho = %.3f, p %s",
  nrow(cor_df), unname(pearson$estimate), format_p(pearson$p.value),
  unname(spearman$estimate), format_p(spearman$p.value)
)
p_cor <- ggplot(cor_df, aes(x = ECM_Fib, y = Activated_Fib)) +
  geom_point(size = 1.45, alpha = 0.52, color = "#2878B5") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9,
              color = "#C82423", fill = "#F3B6B2", alpha = 0.28) +
  annotate("label", x = Inf, y = Inf, label = cor_label,
           hjust = 1.05, vjust = 1.1, size = 4, linewidth = 0.25) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  labs(
    title = "PD ECM_Fib vs Activated_Fib across spatial spots",
    subtitle = "RCTD weights_full; each point is one Visium spot",
    x = "ECM_Fib proportion",
    y = "Activated_Fib proportion",
    caption = "Spot-level correlation; not adjusted for spatial autocorrelation."
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "#555555"),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(size = 9, color = "#666666", hjust = 0)
  )
ggsave("PD_ECM_Fib_Activated_Fib_spot_correlation.pdf", p_cor, width = 6.7, height = 5.3)

cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
