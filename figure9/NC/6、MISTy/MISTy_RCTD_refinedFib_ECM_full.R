#!/usr/bin/env Rscript

options(width = 220)
.libPaths(c("G:/gurobi/Rlib45", .libPaths()))

suppressPackageStartupMessages({
  library(mistyR)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

set.seed(20260823)

spatial_file <- "stRNA_RCTD_refinedFib_ECM.rds"
weights_file <- "RCTD_refinedFib_ECM_weights_full.csv"
results_dir <- "MISTy_RCTD_refinedFib_ECM_full"
log_file <- "MISTy_RCTD_refinedFib_ECM_full.log"

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
write.csv(name_map, "MISTy_celltype_name_map_refinedFib_ECM.csv", row.names = FALSE)

input_summary <- data.frame(
  celltype = colnames(rctd_mat),
  mean_proportion = colMeans(rctd_mat),
  sd_proportion = apply(rctd_mat, 2, sd),
  min_proportion = apply(rctd_mat, 2, min),
  max_proportion = apply(rctd_mat, 2, max),
  row.names = NULL
)
write.csv(input_summary, "MISTy_RCTD_refinedFib_ECM_input_summary.csv", row.names = FALSE)

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
  seed = 20260823,
  cv.folds = 10,
  cached = FALSE
)
saveRDS(misty_result, "MISTy_RCTD_refinedFib_ECM_run_return.rds")

misty_collect <- collect_results(results_dir)
saveRDS(misty_collect, "MISTy_RCTD_refinedFib_ECM_collect.rds")

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
  "MISTy_RCTD_refinedFib_ECM_improvement.pdf",
  plot_improvement_stats(misty_collect),
  width = 9,
  height = 6
)

view_names <- unique(misty_collect$importances.aggregated$view)
for (view_name in view_names) {
  safe_view <- gsub("[^A-Za-z0-9_.-]", "_", view_name)
  tryCatch(
    save_misty_plot(
      paste0("MISTy_RCTD_refinedFib_ECM_interactions_", safe_view, ".pdf"),
      plot_interaction_heatmap(misty_collect, view = view_name, cutoff = 0),
      width = 9,
      height = 8
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
ggsave("ECM_Fib_Activated_Fib_spatial_MISTy_input.pdf", p_focus, width = 10, height = 5)

cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
