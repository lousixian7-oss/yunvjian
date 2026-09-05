suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
})

project_dir <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final"
object_file <- file.path(project_dir,
  "results/objects/fibroblast_revision_final_no_unassigned.rds")
out_file <- file.path(project_dir,
  "results/tables/figure6_revised_GK/JK_C3_AGT_state_pseudobulk_edgeR.csv")

fibro <- readRDS(object_file)
DefaultAssay(fibro) <- "RNA"
fibro[["RNA"]] <- JoinLayers(fibro[["RNA"]])
counts <- GetAssayData(fibro, assay = "RNA", layer = "counts")
meta <- fibro[[]]
meta$cell <- rownames(meta)
meta$state_final <- as.character(meta$state_final)

states <- c("PI16_Fib", "Activated_Fib", "ECM_Fib",
            "Inflammatory_Fib", "Adventitial_Fib")
targets <- c("C3", "AGT")

run_state <- function(state_name) {
  sm <- meta %>%
    filter(state_final == state_name) %>%
    count(sample, group, dataset, name = "n_cells") %>%
    filter(n_cells >= 10) %>%
    arrange(sample)
  if (nrow(sm) < 6 || any(table(sm$group) < 3)) return(NULL)

  pb <- vapply(sm$sample, function(s) {
    cells <- meta$cell[meta$state_final == state_name & meta$sample == s]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(counts)))
  rownames(pb) <- rownames(counts)
  colnames(pb) <- sm$sample

  sm$group <- factor(sm$group, levels = c("NC", "PD"))
  sm$dataset <- droplevels(factor(sm$dataset))
  if (nlevels(sm$dataset) > 1) {
    design <- model.matrix(~ dataset + group, data = sm)
    formula_used <- "~ dataset + group"
  } else {
    design <- model.matrix(~ group, data = sm)
    formula_used <- "~ group"
  }

  y <- DGEList(counts = pb)
  keep <- filterByExpr(y, design = design)
  keep[match(targets, rownames(y), nomatch = 0)] <- TRUE
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "groupPD")
  full <- topTags(qlf, n = Inf, sort.by = "none")$table

  bind_rows(lapply(targets, function(g) {
    if (!g %in% rownames(full)) return(NULL)
    data.frame(
      state = state_name,
      gene = g,
      logFC_PD_vs_NC = full[g, "logFC"],
      p_value = full[g, "PValue"],
      FDR_all_genes_within_state = full[g, "FDR"],
      n_NC = sum(sm$group == "NC"),
      n_PD = sum(sm$group == "PD"),
      datasets = paste(levels(sm$dataset), collapse = ";"),
      minimum_cells_per_sample_state = 10,
      model = formula_used,
      stringsAsFactors = FALSE
    )
  }))
}

res <- bind_rows(lapply(states, run_state)) %>%
  group_by(gene) %>%
  mutate(FDR_across_states_within_gene = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(FDR_global_targeted_tests = p.adjust(p_value, method = "BH"),
         direction = ifelse(logFC_PD_vs_NC > 0, "PD higher", "PD lower")) %>%
  arrange(gene, FDR_across_states_within_gene)

write.csv(res, out_file, row.names = FALSE)
print(res, n = Inf)
