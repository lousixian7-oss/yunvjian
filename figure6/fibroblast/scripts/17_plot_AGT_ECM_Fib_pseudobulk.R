.libPaths(c("C:/Users/32266/AppData/Local/R/win-library/4.5", .libPaths()))
pkgs <- c("Seurat", "Matrix", "edgeR", "dplyr", "ggplot2")
stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(edgeR); library(dplyr); library(ggplot2)})
root <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final"
out <- file.path(root, "results/figures/AGT_ECM_Fib_pseudobulk_20260831")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
write.csv(data.frame(package = pkgs, version = vapply(pkgs, function(p) as.character(packageVersion(p)), character(1))), file.path(out, "package_versions.csv"), row.names = FALSE)
cat("Loading original fibroblast object...\n")
fib <- readRDS(file.path(root, "results/objects/fibroblast_revision_final_no_unassigned.rds"))
DefaultAssay(fib) <- "RNA"
fib[["RNA"]] <- JoinLayers(fib[["RNA"]])
counts <- GetAssayData(fib, assay = "RNA", layer = "counts")
meta <- fib[[]]
stopifnot(identical(colnames(counts), rownames(meta)))
meta$cell <- rownames(meta)
meta$state_final <- as.character(meta$state_final)
audit <- meta %>% distinct(sample, group, dataset) %>%
  left_join(meta %>% filter(state_final == "ECM_Fib") %>% count(sample, name = "n_ECM_Fib"), by = "sample") %>%
  mutate(n_ECM_Fib = coalesce(n_ECM_Fib, 0L), included = n_ECM_Fib >= 10) %>% arrange(sample)
stopifnot(!anyDuplicated(audit$sample))
write.csv(audit, file.path(out, "sample_inclusion_audit.csv"), row.names = FALSE)
sm <- audit %>% filter(included)
sm$group <- factor(sm$group, levels = c("NC", "PD"))
sm$dataset <- droplevels(factor(sm$dataset))
cat("Eligible samples:\n"); print(table(sm$group, sm$dataset))
pb <- vapply(sm$sample, function(s) {
  cells <- meta$cell[meta$sample == s & meta$state_final == "ECM_Fib"]
  Matrix::rowSums(counts[, cells, drop = FALSE])
}, numeric(nrow(counts)))
rownames(pb) <- rownames(counts); colnames(pb) <- sm$sample
stopifnot(all(pb >= 0), all(pb == round(pb)))
design <- model.matrix(~ dataset + group, sm)
stopifnot(qr(design)$rank == ncol(design))
y <- DGEList(counts = pb)
keep <- filterByExpr(y, design = design)
keep[match(c("C3", "AGT"), rownames(y), nomatch = 0)] <- TRUE
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
cat("Refitting the original ECM_Fib model...\n")
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "groupPD")
full <- topTags(test, n = Inf, sort.by = "none")$table
agt <- full["AGT", , drop = FALSE]
previous <- read.csv(file.path(root, "results/tables/figure6_revised_GK/JK_C3_AGT_state_pseudobulk_edgeR.csv")) %>% filter(state == "ECM_Fib", gene == "AGT")
stopifnot(nrow(previous) == 1)
check <- data.frame(metric = c("log2FC", "PValue", "FDR_all_genes_within_state"),
  previous = c(previous$logFC_PD_vs_NC, previous$p_value, previous$FDR_all_genes_within_state),
  reproduced = c(agt$logFC, agt$PValue, agt$FDR))
check$abs_difference <- abs(check$previous - check$reproduced)
write.csv(check, file.path(out, "statistics_reproduction_check.csv"), row.names = FALSE)
print(check)
if (any(check$abs_difference > 1e-6)) stop("Reproduced statistics differ; review before plotting.")
sm$AGT_raw_count <- as.numeric(pb["AGT", ])
sm$library_size_filtered <- y$samples$lib.size
sm$TMM_factor <- y$samples$norm.factors
sm$AGT_TMM_CPM <- as.numeric(cpm(y, log = FALSE)["AGT", ])
sm$AGT_log2_CPM <- as.numeric(cpm(y, log = TRUE, prior.count = 2)["AGT", ])
write.csv(sm, file.path(out, "AGT_ECM_Fib_pseudobulk_sample_values.csv"), row.names = FALSE)
write.csv(cbind(gene = "AGT", state = "ECM_Fib", agt, n_NC = sum(sm$group == "NC"), n_PD = sum(sm$group == "PD"), genes_tested = nrow(full)), file.path(out, "AGT_ECM_Fib_pseudobulk_statistics.csv"), row.names = FALSE)
fdr <- agt$FDR
star <- if (fdr < .001) "***" else if (fdr < .01) "**" else if (fdr < .05) "*" else "ns"
yr <- range(sm$AGT_log2_CPM); span <- max(diff(yr), 1)
bracket <- yr[2] + .12 * span
colors <- c(GSE152042 = "#E58B27", GSE164241 = "#00A582", GSE171213 = "#3B7DDD")
caption <- paste0("Inclusion: >=10 ECM_Fib cells per GSM; samples below threshold are not shown.\n",
  "edgeR quasi-likelihood model: ~ Dataset + group; TMM-normalized counts.\n",
  "BH FDR across all ", format(nrow(full), big.mark = ","), " tested genes in ECM_Fib; * FDR < 0.05.")
p <- ggplot(sm, aes(group, AGT_log2_CPM)) +
  geom_boxplot(aes(fill = group, colour = group), width = .55, linewidth = .65, outlier.shape = NA, alpha = .20, show.legend = FALSE) +
  scale_fill_manual(values = c(NC = "#37B9C7", PD = "#FA756E")) +
  scale_colour_manual(values = c(NC = "#37B9C7", PD = "#FA756E", colors), breaks = levels(sm$dataset), name = "Dataset") +
  geom_point(aes(colour = dataset), position = position_jitter(width = .115, height = 0, seed = 831), size = 2.8, alpha = .95) +
  annotate("segment", x = 1, xend = 2, y = bracket, yend = bracket, linewidth = .55) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = bracket, yend = bracket - .035 * span, linewidth = .55) +
  annotate("text", x = 1.5, y = bracket + .043 * span, label = star, size = 6, fontface = "bold") +
  scale_x_discrete(labels = c(NC = paste0("NC\n(n = ", sum(sm$group == "NC"), " GSMs)"), PD = paste0("PD\n(n = ", sum(sm$group == "PD"), " GSMs)"))) +
  scale_y_continuous(expand = expansion(mult = c(.06, .16))) +
  labs(title = "AGT pseudobulk expression in ECM_Fib",
    subtitle = sprintf("Each point represents one GSM sample\nDataset-adjusted FDR = %.4f  |  log2FC (PD vs NC) = %.2f", fdr, agt$logFC),
    x = "Disease group", y = "Sample-level AGT expression (log2 CPM)", caption = caption) +
  theme_classic(base_size = 13, base_family = "Arial") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = .5),
    plot.subtitle = element_text(size = 11, hjust = .5, lineheight = 1.12, margin = margin(b = 5)),
    axis.title = element_text(face = "bold", size = 12), axis.text = element_text(colour = "#222222", size = 11),
    axis.title.x = element_text(margin = margin(t = 7)),
    panel.border = element_rect(fill = NA, colour = "#555555", linewidth = .65),
    axis.line = element_blank(), legend.position = "top", legend.justification = "left",
    legend.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 10),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0, "pt"),
    plot.caption = element_text(size = 9, hjust = 0, lineheight = 1.15, margin = margin(t = 10)),
    plot.margin = margin(10, 12, 9, 10))
stem <- file.path(out, "AGT_ECM_Fib_pseudobulk_sample_level")
ggsave(paste0(stem, ".png"), p, width = 6.2, height = 6.4, dpi = 400, bg = "white")
ggsave(paste0(stem, ".pdf"), p, width = 6.2, height = 6.4, device = cairo_pdf, bg = "white")
ggsave(paste0(stem, ".tiff"), p, width = 6.2, height = 6.4, dpi = 600, compression = "lzw", bg = "white")
writeLines(c("AGT in ECM_Fib: sample-level pseudobulk plot",
  "Raw RNA counts were summed within each GSM x ECM_Fib combination. Only combinations with >=10 cells were retained.",
  "The original script's filterByExpr procedure was preserved, with C3 and AGT retained, followed by TMM normalization.",
  "Boxes: median and interquartile range; whiskers: most extreme observations within 1.5 IQR. All eligible samples are plotted.",
  "Display: edgeR log2 CPM with prior.count=2, not residualized for Dataset. Inference: raw-count edgeR robust QL model ~ dataset + group.",
  "FDR is BH-adjusted across all tested genes within ECM_Fib, not just AGT or a small set of targeted tests.",
  "The GSM is the aggregation unit; independent-donor equivalence was not verified. No donor random effect is included.",
  paste("Prior-result reproduction: maximum absolute difference", max(check$abs_difference))), file.path(out, "plot_methods_and_interpretation.txt"))
cat("DONE: ", out, "\n")
