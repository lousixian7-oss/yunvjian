args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1L) args[[1]] else "G:/1Yunvjian/0a26.7.7singlecell/8.10"
output_dir <- if (length(args) >= 2L) args[[2]] else file.path(root_dir, "fibroblast_revision_final")
input_rds <- file.path(output_dir, "results/objects/fibroblast_clean_harmony_fixed_programs_preliminary.rds")

figure_dir <- file.path(output_dir, "results/figures")
main_dir <- file.path(figure_dir, "main_figure6")
supp_dir <- file.path(figure_dir, "supplement")
table_dir <- file.path(output_dir, "results/tables")
object_dir <- file.path(output_dir, "results/objects")
log_dir <- file.path(output_dir, "logs")
invisible(lapply(c(main_dir, supp_dir, table_dir, object_dir, log_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

log_con <- file(file.path(log_dir, "02_sample_level_mechanism_and_sensitivity.log"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_con)
}, add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(limma)
  library(slingshot)
  library(SingleCellExperiment)
})

seed <- 20260821L
set.seed(seed)

save_plot <- function(path, plot, width, height) {
  ggsave(path, plot, width = width, height = height, device = cairo_pdf, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", path), plot, width = width, height = height,
         dpi = 300, limitsize = FALSE)
}

add_score <- function(object, genes, score_name) {
  genes <- intersect(genes, rownames(object))
  if (length(genes) < 2L) {
    object[[score_name]] <- NA_real_
    return(object)
  }
  object <- AddModuleScore(object, features = list(genes), name = score_name,
                           seed = seed, search = FALSE)
  object[[score_name]] <- object[[paste0(score_name, "1"), drop = TRUE]]
  object
}

cramers_v <- function(tab) {
  test <- suppressWarnings(chisq.test(tab, correct = FALSE))
  sqrt(as.numeric(test$statistic) /
         (sum(tab) * min(nrow(tab) - 1L, ncol(tab) - 1L)))
}

complete_composition <- function(metadata, state_col) {
  sample_meta <- metadata %>% distinct(sample, group, dataset, cohort)
  states <- sort(unique(as.character(metadata[[state_col]])))
  counts <- metadata %>%
    mutate(state = as.character(.data[[state_col]])) %>%
    dplyr::count(sample, state, name = "n")
  tidyr::expand_grid(sample = sample_meta$sample, state = states) %>%
    left_join(counts, by = c("sample", "state")) %>%
    mutate(n = replace_na(n, 0L)) %>%
    left_join(sample_meta, by = "sample") %>%
    group_by(sample) %>%
    mutate(total_fibro = sum(n), fraction = n / total_fibro) %>%
    ungroup()
}

run_composition_limma <- function(comp, minimum_cells, subset_label, formula_string) {
  part <- comp %>% filter(total_fibro >= minimum_cells)
  samples <- part %>%
    distinct(sample, group, dataset, cohort, total_fibro) %>%
    arrange(sample) %>%
    mutate(across(where(is.factor), droplevels))
  count_matrix <- xtabs(n ~ state + sample, data = part)
  count_matrix <- count_matrix[, samples$sample, drop = FALSE]
  prop_matrix <- sweep(count_matrix, 2, colSums(count_matrix), "/")
  transformed <- asin(sqrt(prop_matrix))
  design <- model.matrix(as.formula(formula_string), samples)
  fit <- eBayes(lmFit(transformed, design), robust = TRUE, trend = TRUE)
  tt <- topTable(fit, coef = "groupPD", number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("state") %>%
    dplyr::rename(effect_arcsin_sqrt = logFC, p_value = P.Value, FDR = adj.P.Val)
  raw_summary <- part %>%
    group_by(state, group) %>%
    summarise(mean_fraction = mean(fraction), median_fraction = median(fraction),
              n_samples = n_distinct(sample), .groups = "drop") %>%
    pivot_wider(names_from = group,
                values_from = c(mean_fraction, median_fraction, n_samples))
  tt %>%
    select(state, effect_arcsin_sqrt, t, p_value, FDR, B) %>%
    left_join(raw_summary, by = "state") %>%
    mutate(subset = subset_label, minimum_cells = minimum_cells,
           formula = formula_string, n_samples_total = nrow(samples))
}

run_program_limma <- function(sample_program, programs, formula_string, subset_label) {
  meta <- sample_program %>%
    distinct(sample, group, dataset, cohort, total_fibro) %>%
    arrange(sample) %>%
    mutate(across(where(is.factor), droplevels))
  matrix <- t(as.matrix(sample_program[match(meta$sample, sample_program$sample), programs, drop = FALSE]))
  colnames(matrix) <- meta$sample
  design <- model.matrix(as.formula(formula_string), meta)
  fit <- eBayes(lmFit(matrix, design), robust = TRUE, trend = TRUE)
  topTable(fit, coef = "groupPD", number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("program") %>%
    dplyr::rename(effect = logFC, p_value = P.Value, FDR = adj.P.Val) %>%
    select(program, effect, t, p_value, FDR, B) %>%
    mutate(subset = subset_label, formula = formula_string,
           n_samples_total = nrow(meta))
}

message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
fibro <- readRDS(input_rds)
DefaultAssay(fibro) <- "RNA"
fibro[["RNA"]] <- JoinLayers(fibro[["RNA"]])

annotation <- data.frame(
  cluster_final = as.character(0:8),
  final_state = c(
    "Activated_Fib", "Inflammatory_Fib", "ECM_Fib", "ECM_Fib",
    "Activated_Fib", "Adventitial_Fib", "PI16_Fib",
    "Inflammatory_Fib", "Unassigned_Fib"
  ),
  evidence = c(
    "Activated module with CXCL6/PTGS2 and stress-high flag; retained as activated-program state",
    "CXCL13/CXCL2/IFITM1 inflammatory chemokine program",
    "IGFBP2/COL18A1 matrix-supporting fibroblast state",
    "ECM and activated modules; MMP13/COL1 program; dataset-sensitive",
    "POSTN/COL11A1/TAGLN/ASPN matrix activation",
    "APOD/PTGDS/CFD/ABCA8",
    "PI16/MFAP5/DPT/CD34",
    "CCL19/CXCL2/IL6 inflammatory chemokine program",
    "Rare F13A1/ITIH5 state with low program margin and sample dominance"
  ),
  confidence = c(
    "Moderate-stress-sensitive", "Moderate", "Low-mixed", "Moderate-dataset-sensitive",
    "Moderate", "Moderate", "High", "Moderate-sample-sensitive", "Unassigned"
  ),
  stringsAsFactors = FALSE
)
write.csv(annotation, file.path(table_dir, "13_final_cluster_annotation_dictionary.csv"), row.names = FALSE)

state_map <- setNames(annotation$final_state, annotation$cluster_final)
confidence_map <- setNames(annotation$confidence, annotation$cluster_final)
fibro$state_final <- factor(
  unname(state_map[as.character(fibro$cluster_final)]),
  levels = c("PI16_Fib", "Activated_Fib", "ECM_Fib", "Adventitial_Fib",
             "Inflammatory_Fib", "Unassigned_Fib")
)
fibro$state_confidence <- unname(confidence_map[as.character(fibro$cluster_final)])

inflammatory_genes <- c("ICAM1", "CCL2", "CXCL2", "CXCL6", "CXCL12", "CXCL13", "IL6")
fibro <- add_score(fibro, inflammatory_genes, "score_Inflammatory_fibroblast")
fibro <- add_score(
  fibro,
  c("POSTN", "TGFB1", "TGFB2", "TGFB3", "COL1A1", "COL3A1", "COL5A1",
    "MMP2", "MMP14", "TIMP1"),
  "score_ECM_TGFB_remodeling"
)
fibro <- add_score(
  fibro,
  c("C3", "COL1A1", "CD44", "POSTN", "TNC", "SERPINE1"),
  "score_C3_CD44_activation"
)

metadata <- fibro[[]] %>% mutate(cell = rownames(fibro[[]]))
state_counts <- as.data.frame(table(metadata$dataset, metadata$state_final))
colnames(state_counts) <- c("dataset", "state", "n")
state_table <- xtabs(n ~ dataset + state, state_counts)
write.csv(as.data.frame.matrix(state_table),
          file.path(table_dir, "14_dataset_by_final_state_contingency.csv"))
write.csv(
  data.frame(association = "dataset x final fibroblast state",
             cramers_v = cramers_v(state_table)),
  file.path(table_dir, "15_final_state_dataset_cramers_v.csv"), row.names = FALSE
)

message("Computing sample-level composition and program scores...")
composition <- complete_composition(metadata, "state_final")
write.csv(composition, file.path(table_dir, "16_sample_level_fibroblast_composition.csv"),
          row.names = FALSE)

program_cols <- c(
  "score_PI16_fibroblast", "score_Activated_fibroblast", "score_ECM_fibroblast",
  "score_Adventitial_fibroblast", "score_Myofibroblast",
  "score_Inflammatory_fibroblast", "score_Stress", "score_Contamination",
  "score_ECM_TGFB_remodeling", "score_C3_CD44_activation"
)
sample_program <- metadata %>%
  group_by(sample, group, dataset, cohort) %>%
  summarise(total_fibro = n(), across(all_of(program_cols), mean), .groups = "drop")
write.csv(sample_program, file.path(table_dir, "17_sample_level_module_scores.csv"),
          row.names = FALSE)

composition_stats <- bind_rows(
  run_composition_limma(composition, 50L,
                        "Integrated eligible samples", "~ dataset + group"),
  run_composition_limma(composition, 100L,
                        "Integrated eligible samples sensitivity", "~ dataset + group"),
  run_composition_limma(composition %>% filter(dataset == "GSE164241"), 50L,
                        "GSE164241 cohort-controlled sensitivity", "~ cohort + group"),
  run_composition_limma(composition %>% filter(dataset == "GSE164241"), 100L,
                        "GSE164241 cohort-controlled sensitivity", "~ cohort + group")
)
write.csv(composition_stats, file.path(table_dir, "18_sample_level_composition_statistics.csv"),
          row.names = FALSE)

program_stats <- bind_rows(
  run_program_limma(
    sample_program %>% filter(total_fibro >= 50), program_cols,
    "~ dataset + group", "Integrated eligible samples"
  ),
  run_program_limma(
    sample_program %>% filter(total_fibro >= 100), program_cols,
    "~ dataset + group", "Integrated eligible samples sensitivity"
  ),
  run_program_limma(
    sample_program %>% filter(dataset == "GSE164241", total_fibro >= 50), program_cols,
    "~ cohort + group", "GSE164241 cohort-controlled sensitivity"
  )
)
write.csv(program_stats, file.path(table_dir, "19_sample_level_module_score_statistics.csv"),
          row.names = FALSE)

eligible_counts <- sample_program %>%
  mutate(eligible_50 = total_fibro >= 50, eligible_100 = total_fibro >= 100) %>%
  dplyr::count(dataset, group, eligible_50, eligible_100, name = "n_samples")
write.csv(eligible_counts, file.path(table_dir, "20_eligible_sample_counts.csv"), row.names = FALSE)

message("Building sample-level Figure 6 panels...")
state_colors <- c(
  PI16_Fib = "#D4A017", Activated_Fib = "#D55E00", ECM_Fib = "#7B61A8",
  Adventitial_Fib = "#009E73", Inflammatory_Fib = "#56B4E9", Unassigned_Fib = "#999999"
)
group_colors <- c(NC = "#4DBBD5", PD = "#E64B35")

sample_order <- composition %>%
  distinct(sample, dataset, cohort, group) %>%
  arrange(dataset, cohort, group, sample) %>% pull(sample)
composition_plot_data <- composition %>%
  mutate(sample = factor(sample, levels = sample_order))
p_comp <- ggplot(composition_plot_data,
                 aes(x = sample, y = fraction, fill = state)) +
  geom_col(width = 0.9, color = "white", linewidth = 0.1) +
  facet_grid(~ dataset + group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Fibroblast-state composition by biological sample",
       subtitle = "Each bar represents one GSM sample; denominator is all qualified fibroblasts in that sample",
       x = NULL, y = "Within-sample fraction", fill = "Fibroblast state") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right")
save_plot(file.path(main_dir, "F6C_sample_level_fibroblast_composition.pdf"), p_comp, 17, 6.5)

eligible_program <- sample_program %>% filter(total_fibro >= 50)
p_act <- ggplot(eligible_program,
                aes(x = group, y = score_Activated_fibroblast, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.15) +
  geom_point(aes(shape = dataset), position = position_jitter(width = 0.12),
             size = 2.4, alpha = 0.9) +
  scale_color_manual(values = group_colors) +
  labs(title = "PD-associated enrichment of the activated fibroblast program",
       subtitle = "Sample-level mean module score; samples with at least 50 qualified fibroblasts",
       x = NULL, y = "Activated fibroblast module score", color = NULL, shape = "Dataset") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
save_plot(file.path(main_dir, "F6D_sample_level_activated_program.pdf"), p_act, 7.5, 6)

p_act_dataset <- ggplot(eligible_program,
                        aes(x = group, y = score_Activated_fibroblast, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.15) +
  geom_point(position = position_jitter(width = 0.12), size = 2.2) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_color_manual(values = group_colors) +
  labs(title = "Activated fibroblast program stratified by dataset",
       x = NULL, y = "Sample-level module score", color = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
save_plot(file.path(supp_dir, "S04_activated_program_by_dataset.pdf"), p_act_dataset, 10, 5)

p_gse164 <- ggplot(
  eligible_program %>% filter(dataset == "GSE164241"),
  aes(x = group, y = score_Activated_fibroblast, color = group, shape = cohort)
) +
  geom_boxplot(aes(group = group), outlier.shape = NA, width = 0.55, alpha = 0.15) +
  geom_point(position = position_jitter(width = 0.12), size = 2.5) +
  scale_color_manual(values = group_colors) +
  labs(title = "GSE164241 cohort-controlled sensitivity analysis",
       subtitle = "GSM500 and GSM517 cohorts are distinguished by point shape",
       x = NULL, y = "Activated fibroblast module score", color = NULL, shape = "Cohort") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
save_plot(file.path(supp_dir, "S05_GSE164241_activated_program_cohort_sensitivity.pdf"),
          p_gse164, 8, 6)

heat_data <- state_counts %>%
  group_by(state) %>% mutate(fraction_within_state = n / sum(n)) %>% ungroup()
p_heat <- ggplot(heat_data, aes(x = dataset, y = state, fill = fraction_within_state)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(fraction_within_state, accuracy = 0.1)), size = 3) +
  scale_fill_viridis_c(limits = c(0, 1), labels = scales::percent) +
  labs(title = "Dataset contribution to final fibroblast states",
       subtitle = paste0("Cramer's V = ", sprintf("%.3f", cramers_v(state_table))),
       x = NULL, y = NULL, fill = "Within-state\nfraction") +
  theme_classic(base_size = 10)
save_plot(file.path(supp_dir, "S06_final_state_dataset_association_heatmap.pdf"), p_heat, 8, 5.5)

message("Analyzing AGT/C3 mechanism axis...")
expression_genes <- intersect(c(
  "AGT", "C3", "POSTN", "TGFB1", "TGFB2", "TGFB3", "COL1A1", "COL3A1", "COL5A1",
  "TNC", "SERPINE1", "FAP", "CD44", "MMP2", "MMP14", "TIMP1"
), rownames(fibro))
expression_data <- FetchData(fibro, vars = expression_genes, layer = "data")
metadata[, expression_genes] <- expression_data[rownames(metadata), expression_genes, drop = FALSE]
fibro$AGT_positive <- metadata$AGT > 0
fibro$C3_positive <- metadata$C3 > 0

pi16_cells <- colnames(fibro)[fibro$state_final == "PI16_Fib"]
activated_cells <- colnames(fibro)[fibro$state_final == "Activated_Fib"]

p_agt_feature <- FeaturePlot(fibro, features = "AGT", reduction = "umap_fibro_final",
                             order = TRUE, raster = TRUE) +
  ggtitle("AGT expression in integrated fibroblasts") + theme_classic(base_size = 10)
p_c3_feature <- FeaturePlot(fibro, features = "C3", reduction = "umap_fibro_final",
                            order = TRUE, raster = TRUE) +
  ggtitle("C3 expression in integrated fibroblasts") + theme_classic(base_size = 10)
save_plot(file.path(main_dir, "F6E_F_AGT_C3_featureplots.pdf"),
          p_agt_feature + p_c3_feature, 13, 5.5)

pi16 <- subset(fibro, cells = pi16_cells)
activated <- subset(fibro, cells = activated_cells)
Idents(pi16) <- pi16$group
Idents(activated) <- activated$group
p_agt_vln <- VlnPlot(pi16, features = "AGT", group.by = "group",
                     cols = group_colors, pt.size = 0) +
  ggtitle("AGT in PI16 fibroblasts") + theme_classic(base_size = 11)
p_c3_vln <- VlnPlot(activated, features = "C3", group.by = "group",
                    cols = group_colors, pt.size = 0) +
  ggtitle("C3 in activated fibroblasts") + theme_classic(base_size = 11)
save_plot(file.path(main_dir, "F6E_F_AGT_PI16_C3_Activated_violins.pdf"),
          p_agt_vln + p_c3_vln, 10, 5.5)

mechanism_sample <- metadata %>%
  mutate(
    state_final = as.character(state_final),
    AGT_positive = AGT > 0,
    C3_positive = C3 > 0
  ) %>%
  group_by(sample, group, dataset, cohort, state_final) %>%
  summarise(
    n_cells = n(),
    AGT_mean = mean(AGT), AGT_positive_fraction = mean(AGT_positive),
    C3_mean = mean(C3), C3_positive_fraction = mean(C3_positive),
    ECM_TGFB_score = mean(score_ECM_TGFB_remodeling),
    C3_CD44_score = mean(score_C3_CD44_activation),
    .groups = "drop"
  )
write.csv(mechanism_sample, file.path(table_dir, "21_sample_level_AGT_C3_mechanism_summary.csv"),
          row.names = FALSE)

pi16$AGT_status <- factor(ifelse(pi16$AGT_positive, "AGT+", "AGT-"), levels = c("AGT-", "AGT+"))
activated$C3_status <- factor(ifelse(activated$C3_positive, "C3+", "C3-"), levels = c("C3-", "C3+"))
agt_mechanism_genes <- intersect(c("AGT", "POSTN", "TGFB1", "TGFB2", "TGFB3",
                                   "COL1A1", "COL3A1", "COL5A1", "MMP2", "MMP14", "TIMP1"),
                                 rownames(pi16))
c3_mechanism_genes <- intersect(c("C3", "COL1A1", "CD44", "POSTN", "TNC", "SERPINE1"),
                                rownames(activated))
p_agt_dot <- DotPlot(pi16, features = agt_mechanism_genes, group.by = "AGT_status") +
  RotatedAxis() + labs(title = "ECM/TGF-beta program in AGT+ PI16 fibroblasts", x = NULL, y = NULL) +
  theme_classic(base_size = 9)
p_c3_dot <- DotPlot(activated, features = c3_mechanism_genes, group.by = "C3_status") +
  RotatedAxis() + labs(title = "COL1A1/CD44 program in C3+ activated fibroblasts", x = NULL, y = NULL) +
  theme_classic(base_size = 9)
save_plot(file.path(main_dir, "F6_mechanism_AGT_C3_dotplots.pdf"), p_agt_dot + p_c3_dot, 14, 5.5)

mechanism_cell_summary <- bind_rows(
  pi16[[]] %>%
    group_by(AGT_status) %>%
    summarise(n_cells = n(), n_samples = n_distinct(sample),
              mean_ECM_TGFB_score = mean(score_ECM_TGFB_remodeling), .groups = "drop") %>%
    mutate(analysis = "PI16 fibroblasts: AGT status"),
  activated[[]] %>%
    group_by(C3_status) %>%
    summarise(n_cells = n(), n_samples = n_distinct(sample),
              mean_C3_CD44_score = mean(score_C3_CD44_activation), .groups = "drop") %>%
    mutate(analysis = "Activated fibroblasts: C3 status")
)
write.csv(mechanism_cell_summary, file.path(table_dir, "22_AGT_C3_status_descriptive_summary.csv"),
          row.names = FALSE)

message("Testing PI16-rooted transcriptomic continuity with Slingshot...")
sce <- as.SingleCellExperiment(fibro)
reducedDim(sce, "HARMONY_FINAL") <- Embeddings(fibro, "harmony_fibro_final")[, 1:20, drop = FALSE]
colData(sce)$cluster_final <- fibro$cluster_final
sds <- slingshot(
  sce,
  clusterLabels = "cluster_final",
  reducedDim = "HARMONY_FINAL",
  start.clus = "6",
  allow.breaks = TRUE,
  shrink = TRUE
)
lineages <- slingLineages(sds)
lineage_table <- bind_rows(lapply(seq_along(lineages), function(i) {
  data.frame(lineage = paste0("Lineage", i), order = seq_along(lineages[[i]]),
             cluster = lineages[[i]])
}))
write.csv(lineage_table, file.path(table_dir, "23_PI16_rooted_slingshot_lineages.csv"), row.names = FALSE)

target_lineages <- which(vapply(lineages, function(x) any(x %in% c("0", "4")), logical(1)))
if (length(target_lineages) > 0L) {
  selected_lineage <- target_lineages[[1]]
  pseudotime <- slingPseudotime(sds)[, selected_lineage]
  weights <- slingCurveWeights(sds)[, selected_lineage]
  fibro$PI16_to_Activated_pseudotime <- pseudotime[colnames(fibro)]
  fibro$PI16_to_Activated_lineage_weight <- weights[colnames(fibro)]
  valid <- !is.na(fibro$PI16_to_Activated_pseudotime) &
    fibro$PI16_to_Activated_lineage_weight > 0.5
  trajectory_data <- FetchData(
    fibro,
    vars = c("PI16_to_Activated_pseudotime", "score_PI16_fibroblast",
             "score_Activated_fibroblast", "AGT", "C3", "dataset", "state_final")
  ) %>% filter(valid[rownames(.)])
  correlations <- data.frame(
    variable = c("PI16 module", "Activated module", "AGT", "C3"),
    spearman_rho = c(
      cor(trajectory_data$PI16_to_Activated_pseudotime,
          trajectory_data$score_PI16_fibroblast, method = "spearman", use = "complete.obs"),
      cor(trajectory_data$PI16_to_Activated_pseudotime,
          trajectory_data$score_Activated_fibroblast, method = "spearman", use = "complete.obs"),
      cor(trajectory_data$PI16_to_Activated_pseudotime,
          trajectory_data$AGT, method = "spearman", use = "complete.obs"),
      cor(trajectory_data$PI16_to_Activated_pseudotime,
          trajectory_data$C3, method = "spearman", use = "complete.obs")
    ),
    lineage = paste(lineages[[selected_lineage]], collapse = " -> "),
    n_cells_weight_gt_0.5 = nrow(trajectory_data)
  )
  write.csv(correlations, file.path(table_dir, "24_PI16_to_Activated_pseudotime_correlations.csv"),
            row.names = FALSE)

  p_pt <- FeaturePlot(fibro, features = "PI16_to_Activated_pseudotime",
                      reduction = "umap_fibro_final", order = TRUE, raster = TRUE) +
    ggtitle("PI16-rooted Slingshot pseudotime hypothesis") + theme_classic(base_size = 10)
  save_plot(file.path(supp_dir, "S07_PI16_rooted_slingshot_pseudotime.pdf"), p_pt, 7, 5.5)

  trajectory_deciles <- trajectory_data %>%
    mutate(decile = ntile(PI16_to_Activated_pseudotime, 10)) %>%
    dplyr::count(decile, dataset, name = "n") %>%
    group_by(decile) %>% mutate(fraction = n / sum(n)) %>% ungroup()
  write.csv(trajectory_deciles, file.path(table_dir, "25_trajectory_dataset_composition_by_decile.csv"),
            row.names = FALSE)
} else {
  warning("No PI16-rooted Slingshot lineage reached activated clusters 0 or 4.")
}

saveRDS(fibro, file.path(object_dir, "fibroblast_revision_final.rds"), compress = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_final_analysis.txt"))
message("Final state counts:")
print(table(fibro$state_final))
message("Sample-level activated-program results:")
print(program_stats %>% filter(program == "score_Activated_fibroblast"))
message("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
