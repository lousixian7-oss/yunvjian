options(width = 220)
.libPaths(c("G:/gurobi/Rlib45", .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

set.seed(20260823)

reference_file <- "yjsl_RCTD_fullannotation_refinedFib.rds"
spatial_file <- "spatial_NC_input.rda"
result_file <- "RCTD_refinedFib_ECM_results.rds"
preprocessed_file <- "RCTD_refinedFib_ECM_preprocessed.rds"
spatial_output_file <- "stRNA_RCTD_refinedFib_ECM.rds"
weights_file <- "RCTD_refinedFib_ECM_weights.csv"
weights_full_file <- "RCTD_refinedFib_ECM_weights_full.csv"
annotation_counts_file <- "RCTD_refinedFib_ECM_reference_cell_counts.csv"
run_log <- "RCTD_refinedFib_ECM_run.log"

if (!file.exists(reference_file)) stop("Missing reference file: ", reference_file)
if (!file.exists(spatial_file)) stop("Missing spatial file: ", spatial_file)
if (file.exists(result_file) && Sys.getenv("OVERWRITE_RCTD", "0") != "1") {
  stop("Output exists. Set OVERWRITE_RCTD=1 only if an intentional rerun is required: ", result_file)
}

sink(run_log, split = TRUE)
on.exit({
  while (sink.number() > 0L) sink()
}, add = TRUE)

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R:", R.version.string, "\n")
cat("spacexr:", as.character(packageVersion("spacexr")), "\n")
cat("Mechanism focus: ECM_Fib (replaces the previous PI16_Fib focus)\n")

raw_metadata <- function(object) attr(object, "meta.data", exact = TRUE)

extract_counts <- function(object, assay_name) {
  assays <- attr(object, "assays", exact = TRUE)
  if (!assay_name %in% names(assays)) stop("Missing assay: ", assay_name)
  aa <- attributes(assays[[assay_name]])
  counts <- if (!is.null(aa$layers) && "counts" %in% names(aa$layers)) {
    aa$layers[["counts"]]
  } else {
    aa$counts
  }
  if (is.null(counts)) stop("Missing counts matrix in assay: ", assay_name)

  feature_names <- NULL
  cell_names <- NULL
  if (!is.null(aa$features)) feature_names <- attributes(aa$features)$dimnames[[1L]]
  if (!is.null(aa$cells)) cell_names <- attributes(aa$cells)$dimnames[[1L]]
  if (is.null(feature_names) && !is.null(aa$meta.features)) feature_names <- rownames(aa$meta.features)
  if (is.null(cell_names)) cell_names <- rownames(raw_metadata(object))

  if (length(feature_names) != nrow(counts)) stop("Feature-name/count-row mismatch for assay: ", assay_name)
  if (length(cell_names) != ncol(counts)) stop("Cell-name/count-column mismatch for assay: ", assay_name)
  dimnames(counts) <- list(feature_names, cell_names)
  counts
}

assert_integer_sparse_counts <- function(counts, label) {
  if (!inherits(counts, "sparseMatrix")) stop(label, " counts are not sparse")
  xv <- counts@x
  if (anyNA(xv) || any(xv < 0) || any(xv != round(xv))) stop(label, " counts are not non-negative integers")
  invisible(TRUE)
}

cat("\nLoading refined single-cell reference...\n")
yjsl <- readRDS(reference_file)
yjsl_md <- raw_metadata(yjsl)
if (!"fullannotation" %in% colnames(yjsl_md)) stop("Reference lacks yjsl$fullannotation")

sc_counts <- extract_counts(yjsl, "RNA")
assert_integer_sparse_counts(sc_counts, "Single-cell")
if (!identical(colnames(sc_counts), rownames(yjsl_md))) stop("Reference counts and metadata are not aligned")

cell_types_all <- as.character(yjsl_md$fullannotation)
keep_reference <- !is.na(cell_types_all) & nzchar(cell_types_all) & cell_types_all != "Fibroblasts"
sc_counts <- sc_counts[, keep_reference, drop = FALSE]
cell_types <- factor(cell_types_all[keep_reference])
names(cell_types) <- colnames(sc_counts)
sc_nUMI <- Matrix::colSums(sc_counts)

reference_counts <- sort(table(cell_types), decreasing = TRUE)
if (!all(c("ECM_Fib", "Activated_Fib", "Endothelial cells", "T cells") %in% names(reference_counts))) {
  stop("Required reference annotations are absent")
}
if (any(reference_counts < 25L)) stop("At least one reference type has fewer than 25 cells")
utils::write.csv(
  data.frame(cell_type = names(reference_counts), n_cells = as.integer(reference_counts), row.names = NULL),
  annotation_counts_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
print(reference_counts)
cat("Reference cells retained:", ncol(sc_counts), "\n")

# Deterministic cap at 10,000 cells per type, matching the RCTD default.
indices_by_type <- split(seq_along(cell_types), cell_types)
reference_indices <- unlist(lapply(indices_by_type, function(idx) {
  if (length(idx) > 10000L) sample(idx, 10000L) else idx
}), use.names = FALSE)
reference_indices <- sort(reference_indices)
sc_counts <- sc_counts[, reference_indices, drop = FALSE]
cell_types <- droplevels(cell_types[reference_indices])
names(cell_types) <- colnames(sc_counts)
sc_nUMI <- sc_nUMI[reference_indices]

reference_sce <- SummarizedExperiment(
  assays = list(counts = sc_counts),
  colData = S4Vectors::DataFrame(
    cell_type = cell_types,
    nUMI = as.numeric(sc_nUMI),
    row.names = colnames(sc_counts)
  )
)
cat("Reference cells after the 10,000-per-type cap:", ncol(reference_sce), "\n")

rm(yjsl, sc_counts)
gc()

cat("\nLoading spatial Seurat object...\n")
spatial_env <- new.env(parent = emptyenv())
loaded_names <- load(spatial_file, envir = spatial_env)
if (!"stRNA" %in% loaded_names) stop("Spatial RDA does not contain an object named stRNA")
stRNA <- spatial_env$stRNA
rm(spatial_env)

st_md <- raw_metadata(stRNA)
if (!"Region" %in% colnames(st_md)) stop("stRNA metadata lacks Region")
st_counts_all <- extract_counts(stRNA, "Spatial")
assert_integer_sparse_counts(st_counts_all, "Spatial")

images <- attr(stRNA, "images", exact = TRUE)
if (!"slice1" %in% names(images)) stop("stRNA lacks images$slice1")
image_attrs <- attributes(images[["slice1"]])
centroids <- image_attrs$boundaries[["centroids"]]
centroid_attrs <- attributes(centroids)
coords <- as.data.frame(centroid_attrs$coords)
coord_cells <- centroid_attrs$cells
if (nrow(coords) != length(coord_cells)) stop("Coordinate/cell count mismatch")
rownames(coords) <- coord_cells
colnames(coords)[1:2] <- c("x", "y")
coords <- coords[, c("x", "y"), drop = FALSE]

common_spots <- intersect(colnames(st_counts_all), rownames(coords))
if (length(common_spots) < 100L) stop("Too few spatial spots with coordinates: ", length(common_spots))
st_counts <- st_counts_all[, common_spots, drop = FALSE]
coords <- coords[common_spots, , drop = FALSE]
stopifnot(identical(colnames(st_counts), rownames(coords)))
st_nUMI <- Matrix::colSums(st_counts)
cat("Tissue spots retained:", length(common_spots), "\n")

spatial_spe <- SpatialExperiment(
  assays = list(counts = st_counts),
  colData = S4Vectors::DataFrame(
    nUMI = as.numeric(st_nUMI),
    row.names = colnames(st_counts)
  ),
  spatialCoords = as.matrix(coords)
)

cat("\nPreprocessing RCTD inputs...\n")
rctd_data <- createRctd(
  spatial_experiment = spatial_spe,
  reference_experiment = reference_sce,
  cell_type_col = "cell_type",
  require_int = TRUE,
  ref_n_cells_min = 25,
  ref_n_cells_max = 10000
)
saveRDS(rctd_data, preprocessed_file, compress = "gzip")

cat("\nRunning official RCTD (doublet mode, 12 cores)...\n")
results_spe <- runRctd(
  rctd_data,
  rctd_mode = "doublet",
  max_cores = 12
)
saveRDS(results_spe, result_file, compress = "gzip")

extract_weight_matrix <- function(results, assay_name) {
  if (!assay_name %in% assayNames(results)) stop("RCTD result lacks assay: ", assay_name)
  w <- t(as.matrix(assay(results, assay_name)))
  rs <- rowSums(w)
  if (any(!is.finite(rs)) || any(w < 0)) stop("Invalid RCTD weights in ", assay_name)
  positive <- rs > 0
  w[positive, ] <- w[positive, , drop = FALSE] / rs[positive]
  w[!positive, ] <- 0
  w
}

weights <- extract_weight_matrix(results_spe, "weights")
weights_full <- extract_weight_matrix(results_spe, "weights_full")
if (!all(c("ECM_Fib", "Activated_Fib") %in% colnames(weights))) stop("ECM_Fib or Activated_Fib is absent from RCTD weights")
utils::write.csv(data.frame(spot_id = rownames(weights), weights, check.names = FALSE), weights_file, row.names = FALSE)
utils::write.csv(data.frame(spot_id = rownames(weights_full), weights_full, check.names = FALSE), weights_full_file, row.names = FALSE)

cat("\nAdding RCTD annotations and weights to the spatial Seurat metadata...\n")
result_spots <- rownames(weights)
if (!all(result_spots %in% rownames(st_md))) stop("Some RCTD spots are absent from stRNA metadata")
for (cell_type in colnames(weights)) {
  st_md[[cell_type]] <- NA_real_
  st_md[result_spots, cell_type] <- weights[, cell_type]
}
st_md$RCTD_celltype <- NA_character_
first_type <- as.character(colData(results_spe)$first_type)
names(first_type) <- colnames(results_spe)
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

cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Result spots:", nrow(weights), "\n")
cat("Weight cell types:", paste(colnames(weights), collapse = ", "), "\n")
cat("ECM_Fib summary:\n")
print(summary(weights[, "ECM_Fib"]))
cat("Activated_Fib summary:\n")
print(summary(weights[, "Activated_Fib"]))
