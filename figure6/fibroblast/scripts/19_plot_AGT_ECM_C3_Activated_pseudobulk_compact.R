.libPaths(c("C:/Users/32266/AppData/Local/R/win-library/4.5", .libPaths()))

pkgs <- c("Seurat", "Matrix", "edgeR", "dplyr", "ggplot2", "patchwork")
stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

root <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final"
out <- file.path(root, "results/figures/AGT_ECM_C3_Activated_pseudobulk_20260901")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fib <- readRDS(file.path(root, "results/objects/fibroblast_revision_final_no_unassigned.rds"))
DefaultAssay(fib) <- "RNA"
fib[["RNA"]] <- JoinLayers(fib[["RNA"]])
counts <- GetAssayData(fib, assay = "RNA", layer = "counts")
meta <- fib[[]]
meta$cell <- rownames(meta)
meta$state_final <- as.character(meta$state_final)
stopifnot(identical(colnames(counts), meta$cell))

targeted_file <- file.path(root, "results/tables/figure6_revised_GK/JK_C3_AGT_state_pseudobulk_edgeR.csv")
targeted <- read.csv(targeted_file, check.names = FALSE)
sensitivity_file <- file.path(root, "results/figures/Figure6_AGT_C3_inclusion_sensitivity_20260901/01_AGT_C3_state_pseudobulk_threshold_sensitivity.csv")
sensitivity <- read.csv(sensitivity_file, check.names = FALSE)

run_pb <- function(state_name, gene_name, min_cells = 10L) {
  sm <- meta %>%
    filter(state_final == state_name) %>%
    count(sample, group, dataset, name = "n_cells") %>%
    filter(n_cells >= min_cells) %>%
    arrange(sample)
  stopifnot(sum(sm$group == "NC") >= 3L, sum(sm$group == "PD") >= 3L)

  pb <- vapply(sm$sample, function(s) {
    cells <- meta$cell[meta$state_final == state_name & meta$sample == s]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(counts)))
  rownames(pb) <- rownames(counts)
  colnames(pb) <- sm$sample

  sm$group <- factor(sm$group, levels = c("NC", "PD"))
  sm$dataset <- droplevels(factor(sm$dataset))
  if (nlevels(sm$dataset) > 1L) {
    design_try <- model.matrix(~ dataset + group, data = sm)
  } else {
    design_try <- NULL
  }
  if (!is.null(design_try) && qr(design_try)$rank == ncol(design_try)) {
    design <- design_try
    model_used <- "~ dataset + group"
  } else {
    design <- model.matrix(~ group, data = sm)
    model_used <- "~ group"
  }

  y <- DGEList(counts = pb)
  keep <- filterByExpr(y, design = design)
  keep[match(c("AGT", "C3"), rownames(y), nomatch = 0L)] <- TRUE
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "groupPD")
  full <- topTags(test, n = Inf, sort.by = "none")$table
  stat <- full[gene_name, , drop = FALSE]

  if (min_cells == 10L) {
    reference <- targeted %>% filter(state == state_name, gene == gene_name)
    targeted_global_q <- reference$FDR_global_targeted_tests
    targeted_state_q <- reference$FDR_across_states_within_gene
  } else {
    reference <- sensitivity %>%
      filter(state == state_name, gene == gene_name, .data$min_cells == .env$min_cells)
    targeted_global_q <- NA_real_
    targeted_state_q <- NA_real_
  }
  stopifnot(nrow(reference) == 1L)
  reproduced <- c(stat$logFC, stat$PValue, stat$FDR)
  expected <- c(reference$logFC_PD_vs_NC, reference$p_value, reference$FDR_all_genes_within_state)
  stopifnot(max(abs(reproduced - expected)) < 1e-6)

  sm$raw_count <- as.numeric(pb[gene_name, ])
  sm$log2_CPM <- as.numeric(cpm(y, log = TRUE, prior.count = 2)[gene_name, ])
  list(
    samples = sm,
    stat = data.frame(
      state = state_name, gene = gene_name, minimum_cells = min_cells,
      n_NC = sum(sm$group == "NC"), n_PD = sum(sm$group == "PD"),
      datasets = paste(levels(sm$dataset), collapse = ";"), model = model_used,
      genes_tested = nrow(full), logFC_PD_vs_NC = stat$logFC,
      p_value = stat$PValue, FDR_all_genes_within_state = stat$FDR,
      FDR_global_targeted_tests = targeted_global_q,
      FDR_across_states_within_gene = targeted_state_q
    )
  )
}

agt <- run_pb("ECM_Fib", "AGT", 10L)
c3 <- run_pb("Activated_Fib", "C3", 10L)

write.csv(bind_rows(agt$stat, c3$stat), file.path(out, "pseudobulk_statistics.csv"), row.names = FALSE)
write.csv(bind_rows(mutate(agt$samples, state = "ECM_Fib", gene = "AGT"),
                    mutate(c3$samples, state = "Activated_Fib", gene = "C3")),
          file.path(out, "pseudobulk_sample_values.csv"), row.names = FALSE)

group_cols <- c(NC = "#00BFC4", PD = "#F8766D")
dataset_cols <- c(GSE152042 = "#E58B27", GSE164241 = "#00A582", GSE171213 = "#3B7DDD")

make_plot <- function(x, tag, title, subtitle, ylab, star_basis = c("transcriptome", "targeted")) {
  star_basis <- match.arg(star_basis)
  sm <- x$samples
  st <- x$stat
  q_star <- if (star_basis == "transcriptome") st$FDR_all_genes_within_state else st$FDR_global_targeted_tests
  star <- if (q_star < .001) "***" else if (q_star < .01) "**" else if (q_star < .05) "*" else "ns"
  yr <- range(sm$log2_CPM)
  span <- max(diff(yr), 1)
  bracket <- yr[2] + 0.10 * span

  ggplot(sm, aes(group, log2_CPM)) +
    geom_boxplot(aes(fill = group, colour = group), width = .56, linewidth = .65,
                 outlier.shape = NA, alpha = .20, show.legend = FALSE) +
    scale_fill_manual(values = group_cols) +
    scale_colour_manual(values = c(group_cols, dataset_cols), breaks = levels(sm$dataset),
                        name = "Dataset") +
    geom_point(aes(colour = dataset), position = position_jitter(width = .11, height = 0, seed = 901),
               size = 2.7, alpha = .95) +
    annotate("segment", x = 1, xend = 2, y = bracket, yend = bracket, linewidth = .55) +
    annotate("segment", x = c(1, 2), xend = c(1, 2), y = bracket,
             yend = bracket - .035 * span, linewidth = .55) +
    annotate("text", x = 1.5, y = bracket + .04 * span, label = star,
             size = 5.3, fontface = "bold") +
    scale_x_discrete(labels = c(
      NC = paste0("NC\n(n = ", sum(sm$group == "NC"), " GSMs)"),
      PD = paste0("PD\n(n = ", sum(sm$group == "PD"), " GSMs)")
    )) +
    scale_y_continuous(expand = expansion(mult = c(.06, .17))) +
    labs(tag = tag, title = title, subtitle = subtitle, x = "Disease group", y = ylab) +
    theme_classic(base_size = 11.5, base_family = "Arial") +
    theme(
      plot.tag = element_text(face = "bold", size = 20),
      plot.tag.position = c(.01, .99),
      plot.title = element_text(size = 13.5, face = "bold", hjust = .5),
      plot.subtitle = element_text(size = 9.5, hjust = .5, lineheight = 1.1),
      axis.title = element_text(face = "bold", size = 10.5),
      axis.text = element_text(colour = "#222222", size = 10),
      panel.border = element_rect(fill = NA, colour = "#555555", linewidth = .6),
      axis.line = element_blank(),
      legend.position = "top", legend.justification = "left",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 9),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(0, "pt"),
      plot.margin = margin(8, 9, 7, 8)
    )
}

p_agt <- make_plot(
  agt, "E", "AGT expression in ECM_Fib",
  sprintf("Dataset-adjusted pseudobulk: FDR = %.4f; log2FC = %.2f",
          agt$stat$FDR_all_genes_within_state, agt$stat$logFC_PD_vs_NC),
  "AGT expression (log2 CPM)", "transcriptome"
)

p_c3 <- make_plot(
  c3, "F", "C3 expression in Activated_Fib",
  sprintf("Targeted BH q = %.4f; transcriptome-wide FDR = %.4f; log2FC = %.2f\nEligible GSMs: GSE164241",
          c3$stat$FDR_global_targeted_tests, c3$stat$FDR_all_genes_within_state,
          c3$stat$logFC_PD_vs_NC),
  "C3 expression (log2 CPM)", "targeted"
)

save_plot <- function(p, stem, width = 5.2, height = 5.0) {
  ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height, dpi = 400, bg = "white")
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(out, paste0(stem, ".tiff")), p, width = width, height = height,
         dpi = 600, compression = "lzw", bg = "white")
}

save_plot(p_agt, "Figure6E_AGT_ECM_Fib_pseudobulk_compact")
save_plot(p_c3, "Figure6F_C3_Activated_Fib_pseudobulk_compact")
pair <- p_agt + p_c3 + plot_layout(ncol = 2)
save_plot(pair, "Figure6EF_AGT_ECM_C3_Activated_pseudobulk", width = 10.4, height = 5.0)

# Diagnostic multi-dataset sensitivity version. The >=3 threshold retains sparse
# GSE171213 states but is not substituted for the primary >=10 rule solely to
# improve significance.
agt_mc <- run_pb("ECM_Fib", "AGT", 3L)
c3_mc <- run_pb("Activated_Fib", "C3", 3L)
write.csv(bind_rows(agt_mc$stat, c3_mc$stat),
          file.path(out, "multidataset_min3_sensitivity_statistics.csv"), row.names = FALSE)
write.csv(bind_rows(mutate(agt_mc$samples, state = "ECM_Fib", gene = "AGT"),
                    mutate(c3_mc$samples, state = "Activated_Fib", gene = "C3")),
          file.path(out, "multidataset_min3_sensitivity_sample_values.csv"), row.names = FALSE)

p_agt_mc <- make_plot(
  agt_mc, "E", "AGT expression in ECM_Fib",
  sprintf("Multi-dataset sensitivity (>=3 cells/GSM): FDR = %.4f; log2FC = %.2f",
          agt_mc$stat$FDR_all_genes_within_state, agt_mc$stat$logFC_PD_vs_NC),
  "AGT expression (log2 CPM)", "transcriptome"
)
p_c3_mc <- make_plot(
  c3_mc, "F", "C3 expression in Activated_Fib",
  sprintf("Multi-dataset sensitivity (>=3 cells/GSM): FDR = %.4f; log2FC = %.2f",
          c3_mc$stat$FDR_all_genes_within_state, c3_mc$stat$logFC_PD_vs_NC),
  "C3 expression (log2 CPM)", "transcriptome"
)
pair_mc <- p_agt_mc + p_c3_mc + plot_layout(ncol = 2)
save_plot(pair_mc, "Figure6EF_multidataset_min3_sensitivity", width = 10.4, height = 5.0)

write.csv(data.frame(package = pkgs,
                     version = vapply(pkgs, function(p) as.character(packageVersion(p)), character(1))),
          file.path(out, "package_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(c(
  "Figure 6E-F compact pseudobulk panels.",
  "Both panels use the same objective inclusion rule: at least 10 cells of the indicated fibroblast state per GSM.",
  "Raw RNA counts were summed within GSM x state and analyzed with edgeR robust quasi-likelihood models after TMM normalization.",
  "AGT/ECM_Fib: model ~ dataset + group; the star is based on BH FDR across all genes tested in ECM_Fib.",
  "C3/Activated_Fib: only GSE164241 retained both groups at >=10 cells/GSM; model ~ group.",
  "For C3, the star is based on BH adjustment across the six targeted C3/AGT x eligible-state tests (q = 0.02538).",
  "The C3 transcriptome-wide within-Activated_Fib FDR is 0.07179 and is displayed explicitly; the result should be described as targeted and cohort-dependent.",
  "Threshold sensitivity is reported separately in Figure6_AGT_C3_inclusion_sensitivity_20260901."
  ,"A separate >=3-cell multi-dataset sensitivity figure is exported: AGT remains significant after all-gene FDR correction, whereas C3 does not."
), file.path(out, "plot_methods_and_interpretation.txt"))

cat("DONE:", out, "\n")
