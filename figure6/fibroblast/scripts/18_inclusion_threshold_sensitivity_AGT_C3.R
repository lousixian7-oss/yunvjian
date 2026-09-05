.libPaths(c("C:/Users/32266/AppData/Local/R/win-library/4.5", .libPaths()))

pkgs <- c("Seurat", "Matrix", "edgeR", "dplyr")
stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
})

root <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final"
out <- file.path(root, "results/figures/Figure6_AGT_C3_inclusion_sensitivity_20260901")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

obj_file <- file.path(root, "results/objects/fibroblast_revision_final_no_unassigned.rds")
fib <- readRDS(obj_file)
DefaultAssay(fib) <- "RNA"
fib[["RNA"]] <- JoinLayers(fib[["RNA"]])
counts <- GetAssayData(fib, assay = "RNA", layer = "counts")
meta <- fib[[]]
meta$cell <- rownames(meta)
meta$state_final <- as.character(meta$state_final)
stopifnot(identical(colnames(counts), meta$cell))

thresholds <- c(1L, 2L, 3L, 5L, 10L, 20L, 30L, 50L)
targets <- c("AGT", "C3")
state_gene <- data.frame(
  state = c("ECM_Fib", "Activated_Fib"),
  gene = c("AGT", "C3"),
  stringsAsFactors = FALSE
)

run_one <- function(state_name, gene_name, min_cells) {
  sm <- meta %>%
    filter(state_final == state_name) %>%
    count(sample, group, dataset, name = "n_cells") %>%
    filter(n_cells >= min_cells) %>%
    arrange(sample)

  n_nc <- sum(sm$group == "NC")
  n_pd <- sum(sm$group == "PD")
  if (n_nc < 3L || n_pd < 3L) {
    return(data.frame(
      state = state_name, gene = gene_name, min_cells = min_cells,
      n_NC = n_nc, n_PD = n_pd, datasets = paste(sort(unique(sm$dataset)), collapse = ";"),
      model = NA_character_, genes_tested = NA_integer_, logFC_PD_vs_NC = NA_real_,
      p_value = NA_real_, FDR_all_genes_within_state = NA_real_, status = "insufficient samples"
    ))
  }

  pb <- vapply(sm$sample, function(s) {
    cells <- meta$cell[meta$state_final == state_name & meta$sample == s]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(counts)))
  rownames(pb) <- rownames(counts)
  colnames(pb) <- sm$sample

  sm$group <- factor(sm$group, levels = c("NC", "PD"))
  sm$dataset <- droplevels(factor(sm$dataset))
  if (nlevels(sm$dataset) > 1L) {
    design_full <- model.matrix(~ dataset + group, data = sm)
  } else {
    design_full <- NULL
  }
  if (!is.null(design_full) && qr(design_full)$rank == ncol(design_full)) {
      design <- design_full
      model_used <- "~ dataset + group"
    } else {
    design <- model.matrix(~ group, data = sm)
    model_used <- "~ group"
  }

  y <- DGEList(counts = pb)
  keep <- filterByExpr(y, design = design)
  keep[match(targets, rownames(y), nomatch = 0L)] <- TRUE
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "groupPD")
  full <- topTags(qlf, n = Inf, sort.by = "none")$table

  data.frame(
    state = state_name, gene = gene_name, min_cells = min_cells,
    n_NC = n_nc, n_PD = n_pd, datasets = paste(levels(sm$dataset), collapse = ";"),
    model = model_used, genes_tested = nrow(full),
    logFC_PD_vs_NC = full[gene_name, "logFC"],
    p_value = full[gene_name, "PValue"],
    FDR_all_genes_within_state = full[gene_name, "FDR"], status = "ok"
  )
}

grid <- merge(state_gene, data.frame(min_cells = thresholds), all = TRUE)
res <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  run_one(grid$state[i], grid$gene[i], grid$min_cells[i])
})) %>%
  arrange(gene, state, min_cells) %>%
  group_by(gene) %>%
  mutate(BH_targeted_across_thresholds = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(res, file.path(out, "01_AGT_C3_state_pseudobulk_threshold_sensitivity.csv"), row.names = FALSE)

cell_counts <- meta %>%
  filter(state_final %in% state_gene$state) %>%
  count(sample, group, dataset, state_final, name = "n_cells") %>%
  arrange(state_final, dataset, group, sample)
write.csv(cell_counts, file.path(out, "02_sample_state_cell_counts.csv"), row.names = FALSE)

run_within_dataset <- function(state_name, gene_name, dataset_name, min_cells = 3L) {
  sm <- meta %>%
    filter(state_final == state_name, dataset == dataset_name) %>%
    count(sample, group, dataset, name = "n_cells") %>%
    filter(n_cells >= min_cells) %>%
    arrange(sample)
  n_nc <- sum(sm$group == "NC")
  n_pd <- sum(sm$group == "PD")
  if (n_nc < 3L || n_pd < 3L) {
    return(data.frame(state = state_name, gene = gene_name, dataset = dataset_name,
      min_cells = min_cells, n_NC = n_nc, n_PD = n_pd, logFC_PD_vs_NC = NA_real_,
      p_value = NA_real_, FDR_all_genes_within_dataset_state = NA_real_,
      status = "insufficient samples"))
  }
  pb <- vapply(sm$sample, function(s) {
    cells <- meta$cell[meta$state_final == state_name & meta$sample == s]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(counts)))
  rownames(pb) <- rownames(counts)
  colnames(pb) <- sm$sample
  sm$group <- factor(sm$group, levels = c("NC", "PD"))
  design <- model.matrix(~ group, data = sm)
  y <- DGEList(counts = pb)
  keep <- filterByExpr(y, design = design)
  keep[match(targets, rownames(y), nomatch = 0L)] <- TRUE
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "groupPD")
  full <- topTags(qlf, n = Inf, sort.by = "none")$table
  data.frame(state = state_name, gene = gene_name, dataset = dataset_name,
    min_cells = min_cells, n_NC = n_nc, n_PD = n_pd,
    logFC_PD_vs_NC = full[gene_name, "logFC"], p_value = full[gene_name, "PValue"],
    FDR_all_genes_within_dataset_state = full[gene_name, "FDR"], status = "ok")
}

within_grid <- merge(state_gene, data.frame(dataset = sort(unique(meta$dataset))), all = TRUE)
within <- bind_rows(lapply(seq_len(nrow(within_grid)), function(i) {
  run_within_dataset(within_grid$state[i], within_grid$gene[i], within_grid$dataset[i], 3L)
})) %>% arrange(gene, state, dataset)
write.csv(within, file.path(out, "04_within_dataset_min3_cells.csv"), row.names = FALSE)

versions <- data.frame(
  package = pkgs,
  version = vapply(pkgs, function(p) as.character(packageVersion(p)), character(1))
)
write.csv(versions, file.path(out, "package_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(c(
  "Purpose: sensitivity analysis for objective per-GSM fibroblast-state inclusion thresholds.",
  "All sensitivity thresholds (1, 2, 3, 5, 10, 20, 30, 50 cells per GSM x state) are reported; none are selected based on significance.",
  "Thresholds 1-3 are diagnostic sensitivity analyses for sparse states and are not automatically considered reliable primary inclusion rules.",
  "Raw RNA counts were summed within GSM x state and analyzed by edgeR robust quasi-likelihood models.",
  "Dataset was included when the dataset-plus-group design was full rank; otherwise the model was group-only.",
  "FDR_all_genes_within_state is BH correction across all genes tested in that state/threshold analysis.",
  "BH_targeted_across_thresholds is shown only as a conservative sensitivity-family summary and does not replace the all-gene FDR.",
  "Within-dataset models at the diagnostic >=3-cell threshold are reported separately; datasets with fewer than three GSMs per group are labeled insufficient.",
  "GSM is the statistical unit. Cells are not treated as independent replicates."
), file.path(out, "README.txt"))

print(res, n = Inf)
