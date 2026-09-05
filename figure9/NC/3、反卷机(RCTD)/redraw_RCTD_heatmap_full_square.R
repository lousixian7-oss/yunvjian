#!/usr/bin/env Rscript

options(width = 220)
.libPaths(c("G:/gurobi/Rlib45", .libPaths()))

suppressPackageStartupMessages({
  library(pheatmap)
})

weights_file <- "RCTD_refinedFib_ECM_weights_full.csv"
spatial_file <- "spatial_NC_input.rda"
output_pdf <- "RCTD_refinedFib_ECM_region_heatmap_square.pdf"
output_csv <- "RCTD_refinedFib_ECM_region_mean_proportions_full.csv"

stopifnot(file.exists(weights_file), file.exists(spatial_file))

weights_df <- read.csv(weights_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("spot_id" %in% colnames(weights_df))
rownames(weights_df) <- weights_df$spot_id
weights_df$spot_id <- NULL
weights <- as.matrix(weights_df)
storage.mode(weights) <- "double"

spatial_env <- new.env(parent = emptyenv())
loaded_names <- load(spatial_file, envir = spatial_env)
if (!"stRNA" %in% loaded_names) stop("Spatial RDA does not contain stRNA")
stRNA <- spatial_env$stRNA
st_md <- attr(stRNA, "meta.data", exact = TRUE)
if (!"Region" %in% colnames(st_md)) stop("stRNA metadata lacks Region")

common_spots <- intersect(rownames(weights), rownames(st_md))
if (length(common_spots) == 0) stop("No common spots between weights_full and spatial metadata")
weights <- weights[common_spots, , drop = FALSE]
regions <- as.character(st_md[common_spots, "Region"])
keep <- !is.na(regions) & nzchar(regions)
weights <- weights[keep, , drop = FALSE]
regions <- regions[keep]

region_levels <- unique(regions)
avg_mat <- vapply(
  region_levels,
  function(region) colMeans(weights[regions == region, , drop = FALSE], na.rm = TRUE),
  numeric(ncol(weights))
)
rownames(avg_mat) <- colnames(weights)
colnames(avg_mat) <- region_levels

avg_out <- data.frame(Celltype = rownames(avg_mat), avg_mat, check.names = FALSE)
write.csv(avg_out, output_csv, row.names = FALSE)

row_sd <- apply(avg_mat, 1, sd, na.rm = TRUE)
if (any(!is.finite(row_sd) | row_sd == 0)) {
  bad <- rownames(avg_mat)[!is.finite(row_sd) | row_sd == 0]
  stop("Cannot row-scale invariant cell types in weights_full: ", paste(bad, collapse = ", "))
}

pheatmap(
  avg_mat,
  scale = "row",
  color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
  border_color = "white",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  cellwidth = 24,
  cellheight = 24,
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 90,
  main = "cell-type enrichment",
  filename = output_pdf
)

cat("Saved:", output_pdf, "\n")
cat("Saved:", output_csv, "\n")
cat("pDCs region means:", paste(signif(avg_mat["pDCs", ], 4), collapse = ", "), "\n")
cat("Neutrophils region means:", paste(signif(avg_mat["Neutrophils", ], 4), collapse = ", "), "\n")
