options(width = 220)

paths <- c(
  yjsl = "G:/1Yunvjian/0a26.7.7singlecell/8.10/yjsl_clean_allGSE.rds",
  old_merged = "G:/1Yunvjian/0a26.7.7singlecell/8.10/yjsl_fullannotation_fibro.rds",
  fibro = "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final/results/objects/fibroblast_revision_final_no_unassigned.rds"
)

cat("R:", R.version.string, "\n")
pkg_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else "not installed"
}
cat("Seurat:", pkg_version("Seurat"), "\n")
cat("SeuratObject:", pkg_version("SeuratObject"), "\n")

inspect_object <- function(label, path) {
  cat("\n=====", label, "=====\n")
  cat("path:", path, "\n")
  x <- readRDS(path)
  xa <- attributes(x)
  cat("class:", paste(xa$class, collapse = ", "), "\n")
  md <- xa$meta.data
  assays <- xa$assays
  first_assay_attrs <- if (length(assays)) attributes(assays[[1L]]) else list()
  first_layer <- if (length(first_assay_attrs$layers)) first_assay_attrs$layers[[1L]] else first_assay_attrs$counts
  feature_n <- if (!is.null(first_layer)) {
    layer_dim <- attr(first_layer, "Dim")
    if (length(layer_dim)) layer_dim[[1L]] else nrow(first_layer)
  } else NA_integer_
  cat("dim:", feature_n, "x", nrow(md), "\n")
  cat("assays:", paste(names(assays), collapse = ", "), "\n")
  cat("default assay:", xa$active.assay, "\n")
  cat("object version:", as.character(xa$version), "\n")
  for (assay_name in names(assays)) {
    aa <- attributes(assays[[assay_name]])
    layer_names <- if (!is.null(aa$layers)) names(aa$layers) else intersect(c("counts", "data", "scale.data"), names(aa))
    cat("layers", assay_name, ":", paste(layer_names, collapse = ", "), "\n")
  }
  cat("metadata columns (", ncol(md), "):\n", paste(colnames(md), collapse = " | "), "\n", sep = "")
  candidates <- grep("annot|cell.?type|state|label|cluster|class|lineage|full", colnames(md), ignore.case = TRUE, value = TRUE)
  cat("candidate metadata columns:", paste(candidates, collapse = ", "), "\n")
  for (nm in candidates) {
    vals <- as.character(md[[nm]])
    tab <- sort(table(vals, useNA = "ifany"), decreasing = TRUE)
    cat("\n--", nm, "unique=", length(tab), "--\n")
    print(utils::head(tab, 30))
  }
  invisible(list(cells = rownames(md), metadata = md))
}

y <- inspect_object("yjsl", paths[["yjsl"]])
gc()
o <- inspect_object("old_merged", paths[["old_merged"]])
cat("\nyjsl vs old_merged identical cell names:", identical(y$cells, o$cells), "\n")
if ("fullannotation" %in% colnames(o$metadata)) {
  cat("old fullannotation table:\n")
  print(sort(table(o$metadata$fullannotation, useNA = "ifany"), decreasing = TRUE))
}
rm(o)
gc()
f <- inspect_object("fibro", paths[["fibro"]])
cat("\nfibro cells overlapping yjsl:", sum(f$cells %in% y$cells), "/", length(f$cells), "\n")
cat("fibro cells absent from yjsl:", sum(!f$cells %in% y$cells), "\n")
