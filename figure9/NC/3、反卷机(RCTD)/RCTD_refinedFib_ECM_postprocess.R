options(width = 220)
.libPaths(c("G:/gurobi/Rlib45", .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

result_file <- "RCTD_refinedFib_ECM_results.rds"
spatial_file <- "spatial_NC_input.rda"
spatial_output_file <- "stRNA_RCTD_refinedFib_ECM.rds"
weights_file <- "RCTD_refinedFib_ECM_weights.csv"
weights_full_file <- "RCTD_refinedFib_ECM_weights_full.csv"
postprocess_log <- "RCTD_refinedFib_ECM_postprocess.log"

stopifnot(file.exists(result_file), file.exists(spatial_file))
sink(postprocess_log, split = TRUE)
on.exit({ while (sink.number() > 0L) sink() }, add = TRUE)

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
results_spe <- readRDS(result_file)
cat("Result assays:", paste(assayNames(results_spe), collapse = ", "), "\n")

extract_weight_matrix <- function(results, assay_name) {
  if (!assay_name %in% assayNames(results)) stop("RCTD result lacks assay: ", assay_name)
  w <- t(as.matrix(assay(results, assay_name)))
  rs <- rowSums(w)
  if (any(!is.finite(rs)) || any(w < 0)) stop("Invalid RCTD weights in ", assay_name)
  positive <- rs > 0
  w[positive, ] <- w[positive, , drop = FALSE] / rs[positive]
  w[!positive, ] <- 0
  if (any(positive) && max(abs(rowSums(w[positive, , drop = FALSE]) - 1)) > 1e-10) {
    stop("Normalized positive rows do not sum to one in ", assay_name)
  }
  w
}

weights <- extract_weight_matrix(results_spe, "weights")
weights_full <- extract_weight_matrix(results_spe, "weights_full")
if (!all(c("ECM_Fib", "Activated_Fib") %in% colnames(weights))) {
  stop("ECM_Fib or Activated_Fib is absent from RCTD weights")
}
utils::write.csv(data.frame(spot_id = rownames(weights), weights, check.names = FALSE), weights_file, row.names = FALSE)
utils::write.csv(data.frame(spot_id = rownames(weights_full), weights_full, check.names = FALSE), weights_full_file, row.names = FALSE)

spatial_env <- new.env(parent = emptyenv())
loaded_names <- load(spatial_file, envir = spatial_env)
if (!"stRNA" %in% loaded_names) stop("Spatial RDA does not contain stRNA")
stRNA <- spatial_env$stRNA
rm(spatial_env)
st_md <- attr(stRNA, "meta.data", exact = TRUE)
if (!"Region" %in% colnames(st_md)) stop("stRNA metadata lacks Region")

images <- attr(stRNA, "images", exact = TRUE)
image_attrs <- attributes(images[["slice1"]])
centroid_attrs <- attributes(image_attrs$boundaries[["centroids"]])
coords <- as.data.frame(centroid_attrs$coords)
rownames(coords) <- centroid_attrs$cells
colnames(coords)[1:2] <- c("x", "y")
coords <- coords[, c("x", "y"), drop = FALSE]

result_spots <- rownames(weights)
if (!all(result_spots %in% rownames(st_md))) stop("Some RCTD spots are absent from stRNA metadata")
if (!all(result_spots %in% rownames(coords))) stop("Some RCTD spots are absent from spatial coordinates")

for (cell_type in colnames(weights)) {
  st_md[[cell_type]] <- NA_real_
  st_md[result_spots, cell_type] <- weights[, cell_type]
}
st_md$RCTD_celltype <- NA_character_
first_type <- as.character(colData(results_spe)$first_type)
names(first_type) <- colnames(results_spe)
if (!all(names(first_type) %in% rownames(st_md))) stop("Some predicted spots are absent from stRNA metadata")
st_md[names(first_type), "RCTD_celltype"] <- first_type
attr(stRNA, "meta.data") <- st_md
saveRDS(stRNA, spatial_output_file, compress = "gzip")

plot_df <- data.frame(
  spot_id = result_spots,
  x = coords[result_spots, "x"],
  y = coords[result_spots, "y"],
  Region = st_md[result_spots, "Region"],
  ECM_Fib = weights[result_spots, "ECM_Fib"],
  Activated_Fib = weights[result_spots, "Activated_Fib"],
  row.names = NULL,
  check.names = FALSE
)

spatial_plot <- function(feature, title) {
  ggplot(plot_df, aes(x = x, y = y, color = .data[[feature]])) +
    geom_point(size = 2.8) +
    scale_y_reverse() +
    coord_fixed() +
    scale_color_gradientn(colors = c("#F7FBFF", "#FFFF00", "#D7301F"), trans = "sqrt") +
    theme_void(base_size = 13) +
    labs(title = title, color = "RCTD\nproportion") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_ecm <- spatial_plot("ECM_Fib", "ECM_Fib spatial proportion")
p_activated <- spatial_plot("Activated_Fib", "Activated_Fib spatial proportion")
ggsave("ECM_Fib_spatial_RCTD.pdf", p_ecm, width = 7, height = 6)
ggsave("Activated_Fib_spatial_RCTD.pdf", p_activated, width = 7, height = 6)

region_plot <- function(feature) {
  kw <- kruskal.test(plot_df[[feature]] ~ plot_df$Region)
  ggplot(plot_df, aes(x = Region, y = .data[[feature]], fill = Region)) +
    geom_boxplot(outlier.size = 0.4) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
    labs(
      title = paste(feature, "distribution across regions"),
      subtitle = paste0("Kruskal-Wallis p = ", format.pval(kw$p.value, digits = 3)),
      x = NULL,
      y = paste(feature, "proportion")
    )
}

ggsave("ECM_Fib_region_boxplot_RCTD.pdf", region_plot("ECM_Fib"), width = 7.5, height = 6)
ggsave("Activated_Fib_region_boxplot_RCTD.pdf", region_plot("Activated_Fib"), width = 7.5, height = 6)

weights_df <- data.frame(spot_id = rownames(weights), weights, check.names = FALSE)
composition_long <- weights_df %>%
  pivot_longer(cols = -spot_id, names_to = "Celltype", values_to = "Proportion") %>%
  left_join(data.frame(spot_id = result_spots, Region = st_md[result_spots, "Region"]), by = "spot_id")

avg_celltype <- composition_long %>%
  group_by(Region, Celltype) %>%
  summarise(Mean_proportion = mean(Proportion), .groups = "drop")
utils::write.csv(avg_celltype, "RCTD_refinedFib_ECM_region_mean_proportions.csv", row.names = FALSE, fileEncoding = "UTF-8")

p_stack <- ggplot(avg_celltype, aes(x = Region, y = Mean_proportion, fill = Celltype)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.15) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "RCTD cell-type composition across spatial regions", x = NULL, y = "Mean proportion", fill = "Cell type")
ggsave("RCTD_refinedFib_ECM_region_composition.pdf", p_stack, width = 11, height = 7)

heat_df <- avg_celltype %>%
  group_by(Celltype) %>%
  mutate(z = if (sd(Mean_proportion) > 0) as.numeric(scale(Mean_proportion)) else 0) %>%
  ungroup()
p_heat <- ggplot(heat_df, aes(x = Region, y = Celltype, fill = z)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "RCTD cell-type enrichment", x = NULL, y = NULL, fill = "Row z-score")
ggsave("RCTD_refinedFib_ECM_region_heatmap.pdf", p_heat, width = 8.5, height = 7.5)

p_violin <- ggplot(composition_long, aes(x = Region, y = Proportion, fill = Region)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~Celltype, scales = "free_y", ncol = 5) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
  labs(title = "RCTD cell-type proportions across spatial regions", x = NULL, y = "RCTD proportion")
ggsave("RCTD_refinedFib_ECM_region_violin_all.pdf", p_violin, width = 15, height = 12)

ecm_kw <- kruskal.test(plot_df$ECM_Fib ~ plot_df$Region)
activated_kw <- kruskal.test(plot_df$Activated_Fib ~ plot_df$Region)
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Result spots:", nrow(weights), "\n")
cat("Cell types:", paste(colnames(weights), collapse = ", "), "\n")
positive_weight_rows <- rowSums(weights) > 0
cat("Restricted-weight zero rows (RCTD reject spots):", sum(!positive_weight_rows), "\n")
cat("Max absolute row-sum error among positive rows:", max(abs(rowSums(weights[positive_weight_rows, , drop = FALSE]) - 1)), "\n")
cat("ECM_Fib summary:\n")
print(summary(weights[, "ECM_Fib"]))
cat("ECM_Fib Kruskal-Wallis p:", ecm_kw$p.value, "\n")
cat("Activated_Fib summary:\n")
print(summary(weights[, "Activated_Fib"]))
cat("Activated_Fib Kruskal-Wallis p:", activated_kw$p.value, "\n")
cat("Dominant cell-type counts:\n")
print(sort(table(first_type), decreasing = TRUE))
