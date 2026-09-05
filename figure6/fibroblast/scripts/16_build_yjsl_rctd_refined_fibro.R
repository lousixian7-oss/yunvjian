options(width = 220)

input_yjsl <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/yjsl_clean_allGSE.rds"
input_fibro <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final/results/objects/fibroblast_revision_final_no_unassigned.rds"
output_rds <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/yjsl_RCTD_fullannotation_refinedFib.rds"
output_counts <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/yjsl_RCTD_fullannotation_refinedFib_cell_counts.csv"
output_audit <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/yjsl_RCTD_fullannotation_refinedFib_audit.txt"
temp_rds <- paste0(output_rds, ".tmp")

stopifnot(file.exists(input_yjsl), file.exists(input_fibro))
if (file.exists(output_rds)) stop("Refusing to overwrite existing output: ", output_rds)
if (file.exists(temp_rds)) stop("Temporary output already exists; inspect it before rerunning: ", temp_rds)

raw_metadata <- function(object) attr(object, "meta.data", exact = TRUE)

raw_counts <- function(object, assay_name = "RNA") {
  assays <- attr(object, "assays", exact = TRUE)
  if (!assay_name %in% names(assays)) stop("Missing assay: ", assay_name)
  aa <- attributes(assays[[assay_name]])
  if (!is.null(aa$layers) && "counts" %in% names(aa$layers)) return(aa$layers[["counts"]])
  if (!is.null(aa$counts)) return(aa$counts)
  stop("No RNA counts matrix found")
}

matrix_fingerprint <- function(mat) {
  d <- attr(mat, "Dim", exact = TRUE)
  dn <- attr(mat, "Dimnames", exact = TRUE)
  xv <- attr(mat, "x", exact = TRUE)
  if (is.null(d)) d <- dim(mat)
  if (is.null(dn)) dn <- dimnames(mat)
  if (is.null(dn)) dn <- list(NULL, NULL)
  if (is.null(xv)) xv <- as.numeric(mat)
  first_or_missing <- function(v) if (length(v)) v[[1L]] else "<not stored>"
  last_or_missing <- function(v) if (length(v)) v[[length(v)]] else "<not stored>"
  c(
    features = d[[1L]],
    cells = d[[2L]],
    stored_nonzero = length(xv),
    count_sum = sum(xv),
    first_gene = first_or_missing(dn[[1L]]),
    last_gene = last_or_missing(dn[[1L]]),
    first_cell = first_or_missing(dn[[2L]]),
    last_cell = last_or_missing(dn[[2L]])
  )
}

cat("Loading yjsl master object...\n")
yjsl <- readRDS(input_yjsl)
yjsl_md <- raw_metadata(yjsl)
yjsl_counts <- raw_counts(yjsl)
counts_before <- matrix_fingerprint(yjsl_counts)

cat("Loading revised fibroblast object...\n")
fibro <- readRDS(input_fibro)
fibro_md <- raw_metadata(fibro)

required_yjsl <- c("cellType", "nCount_RNA")
required_fibro <- "state_final"
if (!all(required_yjsl %in% colnames(yjsl_md))) stop("yjsl is missing: ", paste(setdiff(required_yjsl, colnames(yjsl_md)), collapse = ", "))
if (!all(required_fibro %in% colnames(fibro_md))) stop("fibro object is missing: ", paste(setdiff(required_fibro, colnames(fibro_md)), collapse = ", "))
if (anyDuplicated(rownames(yjsl_md))) stop("Duplicated cell names in yjsl metadata")
if (anyDuplicated(rownames(fibro_md))) stop("Duplicated cell names in fibro metadata")

count_cells <- attr(yjsl_counts, "Dimnames", exact = TRUE)[[2L]]
if (is.null(count_cells)) count_cells <- colnames(yjsl_counts)
metadata_cells <- rownames(yjsl_md)
count_dim <- attr(yjsl_counts, "Dim", exact = TRUE)
if (length(count_cells) == 0L && count_dim[[2L]] != length(metadata_cells)) {
  stop("RNA counts column number does not match yjsl metadata row number")
}
if (length(count_cells) > 0L && (length(count_cells) != length(metadata_cells) || any(count_cells != metadata_cells))) {
  stop(
    "RNA counts columns and yjsl metadata rows differ; count-only=",
    length(setdiff(count_cells, metadata_cells)),
    ", metadata-only=",
    length(setdiff(metadata_cells, count_cells))
  )
}

fibro_cells <- rownames(fibro_md)
missing_fibro <- setdiff(fibro_cells, rownames(yjsl_md))
if (length(missing_fibro)) stop(length(missing_fibro), " revised fibro cells are absent from yjsl")

fibro_state <- as.character(fibro_md$state_final)
if (anyNA(fibro_state) || any(!nzchar(fibro_state))) stop("state_final contains missing or blank labels")
expected_states <- c("Activated_Fib", "ECM_Fib", "Inflammatory_Fib", "Adventitial_Fib", "PI16_Fib")
if (!all(expected_states %in% fibro_state)) stop("Missing expected refined states: ", paste(setdiff(expected_states, fibro_state), collapse = ", "))

master_index <- match(fibro_cells, rownames(yjsl_md))
if (anyNA(master_index)) stop("Internal barcode matching failure")
matched_broad_types <- as.character(yjsl_md$cellType[master_index])
if (any(matched_broad_types != "Fibroblasts")) {
  stop("Some revised fibro cells are not Fibroblasts in yjsl$cellType: ", paste(names(sort(table(matched_broad_types), decreasing = TRUE)), collapse = ", "))
}

base_annotation <- as.character(yjsl_md$cellType)
if (anyNA(base_annotation) || any(!nzchar(base_annotation))) stop("yjsl$cellType contains missing or blank labels")
refined_subtype <- rep(NA_character_, nrow(yjsl_md))
refined_subtype[master_index] <- fibro_state
fullannotation <- base_annotation
fullannotation[master_index] <- fibro_state

# Keep both spellings: `fullannotation` is the requested field, while
# `full_annotation` remains compatible with the previous analysis script.
yjsl_md$Fibro_subtype_refined <- refined_subtype
yjsl_md$fullannotation <- factor(fullannotation)
yjsl_md$full_annotation <- yjsl_md$fullannotation
yjsl_md$RCTD_annotation_source <- ifelse(
  seq_len(nrow(yjsl_md)) %in% master_index,
  "fibroblast_revision_final_no_unassigned:state_final",
  "yjsl_clean_allGSE:cellType"
)

unresolved_fibro <- sum(as.character(yjsl_md$fullannotation) == "Fibroblasts")
annotation_table <- sort(table(yjsl_md$fullannotation, useNA = "ifany"), decreasing = TRUE)
if (unresolved_fibro != sum(base_annotation == "Fibroblasts") - length(fibro_cells)) {
  stop("Unexpected number of broad/unresolved fibroblasts")
}
if (!all(c("Endothelial cells", "T cells", "Activated_Fib", "ECM_Fib") %in% names(annotation_table))) {
  stop("Required RCTD annotations are missing")
}
if (anyNA(yjsl_md$fullannotation)) stop("fullannotation contains NA")

attr(yjsl, "meta.data") <- yjsl_md
counts_after <- matrix_fingerprint(raw_counts(yjsl))
if (!identical(counts_before, counts_after)) stop("RNA counts fingerprint changed while adding metadata")

cat("Saving new RCTD-ready yjsl object...\n")
saveRDS(yjsl, temp_rds, compress = "gzip")
if (!file.rename(temp_rds, output_rds)) stop("Could not rename temporary RDS to final output")

count_df <- data.frame(
  fullannotation = names(annotation_table),
  n_cells = as.integer(annotation_table),
  row.names = NULL,
  check.names = FALSE
)
utils::write.csv(count_df, output_counts, row.names = FALSE, fileEncoding = "UTF-8")

audit_lines <- c(
  paste("Created:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("Input yjsl:", input_yjsl),
  paste("Input revised fibro:", input_fibro),
  paste("Output:", output_rds),
  paste("Seurat object version slot:", as.character(attr(yjsl, "version", exact = TRUE))),
  paste("Master cells:", nrow(yjsl_md)),
  paste("Refined fibro cells mapped exactly:", length(fibro_cells)),
  paste("Broad Fibroblasts retained (not present in cleaned revised fibro object):", unresolved_fibro),
  paste("RNA features:", counts_before[["features"]]),
  paste("RNA count cells:", counts_before[["cells"]]),
  paste("RNA stored nonzero entries:", counts_before[["stored_nonzero"]]),
  paste("RNA count sum:", counts_before[["count_sum"]]),
  "Counts fingerprint unchanged after metadata update: TRUE",
  "Annotation source: cellType for non-refined cells; state_final for revised fibro cells",
  "Metadata aliases added: fullannotation, full_annotation",
  "",
  "fullannotation cell counts:",
  paste(names(annotation_table), as.integer(annotation_table), sep = "\t")
)
writeLines(audit_lines, output_audit, useBytes = TRUE)

cat("Verifying saved RDS can be read and has the requested metadata...\n")
rm(fibro, fibro_md, yjsl_counts)
gc()
check <- readRDS(output_rds)
check_md <- raw_metadata(check)
if (!all(c("fullannotation", "full_annotation", "Fibro_subtype_refined") %in% colnames(check_md))) stop("Saved RDS is missing new metadata")
if (!identical(as.character(check_md$fullannotation), as.character(yjsl_md$fullannotation))) stop("Saved fullannotation differs from in-memory metadata")
if (!identical(matrix_fingerprint(raw_counts(check)), counts_before)) stop("Saved RDS counts fingerprint differs from input")

cat("\nSUCCESS\n")
cat("Output:", output_rds, "\n")
cat("Cells:", nrow(check_md), "\n")
cat("Refined fibro cells:", length(fibro_cells), "\n")
cat("Broad Fibroblasts retained:", unresolved_fibro, "\n")
print(annotation_table)
