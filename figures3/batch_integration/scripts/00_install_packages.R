args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript 00_install_packages.R <local_library>")
}

local_library <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 1200
)

required <- c("Seurat", "harmony", "patchwork", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0L) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(
    missing,
    lib = local_library,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    type = "binary"
  )
}

still_missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0L) {
  stop("Package installation incomplete: ", paste(still_missing, collapse = ", "))
}

message("All required packages are available.")

