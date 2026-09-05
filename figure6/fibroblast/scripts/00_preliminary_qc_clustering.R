args <- commandArgs(trailingOnly = TRUE)

root_dir <- if (length(args) >= 1L) args[[1]] else "G:/1Yunvjian/0a26.7.7singlecell/8.10"
output_dir <- if (length(args) >= 2L) args[[2]] else file.path(root_dir, "fibroblast_revision_20260821")
input_rds <- file.path(
  root_dir,
  "batch_integration_revision_20260821/results/objects/yjsl_batch_integration_revision_clustered.rds"
)
old_fibro_rds <- file.path(root_dir, "fibro_clean_for_trajectory.rds")

figure_dir <- file.path(output_dir, "results", "figures")
table_dir <- file.path(output_dir, "results", "tables")
object_dir <- file.path(output_dir, "results", "objects")
log_dir <- file.path(output_dir, "logs")
invisible(lapply(
  c(figure_dir, table_dir, object_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

log_file <- file.path(log_dir, "01_fibroblast_harmony_reanalysis.log")
log_connection <- file(log_file, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_connection)
}, add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(limma)
  library(mclust)
})

seed <- 20260821L
set.seed(seed)
message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("Input: ", input_rds)
message(
  "Versions | R ", getRversion(),
  " | Seurat ", packageVersion("Seurat"),
  " | harmony ", packageVersion("harmony"),
  " | limma ", packageVersion("limma")
)

save_plot <- function(filename, plot, width, height) {
  ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    limitsize = FALSE
  )
}

cohort_from_sample <- function(sample_id) {
  case_when(
    grepl("^GSM460", sample_id) ~ "GSM460_cohort",
    grepl("^GSM500", sample_id) ~ "GSM500_cohort",
    grepl("^GSM517", sample_id) ~ "GSM517_cohort",
    grepl("^GSM522", sample_id) ~ "GSM522_cohort",
    TRUE ~ "Other_cohort"
  )
}

composition_long <- function(metadata, state_col, min_cells = 0L) {
  sample_meta <- metadata %>%
    distinct(sample, group, dataset, cohort)
  states <- sort(unique(as.character(metadata[[state_col]])))
  counts <- metadata %>%
    mutate(state = as.character(.data[[state_col]])) %>%
    dplyr::count(sample, state, name = "n")
  complete <- tidyr::expand_grid(
    sample = sample_meta$sample,
    state = states
  ) %>%
    left_join(counts, by = c("sample", "state")) %>%
    mutate(n = replace_na(n, 0L)) %>%
    left_join(sample_meta, by = "sample") %>%
    group_by(sample) %>%
    mutate(total_fibro = sum(n), proportion = n / total_fibro) %>%
    ungroup() %>%
    filter(total_fibro >= min_cells)
  complete
}

run_propeller_style <- function(comp, subset_name, formula_string) {
  samples <- comp %>%
    distinct(sample, group, dataset, cohort, total_fibro) %>%
    arrange(sample) %>%
    mutate(across(where(is.factor), droplevels))
  state_levels <- sort(unique(comp$state))
  count_matrix <- xtabs(n ~ state + sample, data = comp)
  count_matrix <- count_matrix[state_levels, samples$sample, drop = FALSE]
  prop_matrix <- sweep(count_matrix, 2, colSums(count_matrix), "/")
  transformed <- asin(sqrt(prop_matrix))
  design <- model.matrix(as.formula(formula_string), data = samples)
  if (!"groupPD" %in% colnames(design)) {
    stop("Design does not contain groupPD for ", subset_name)
  }
  fit <- eBayes(lmFit(transformed, design), trend = TRUE, robust = TRUE)
  result <- topTable(fit, coef = "groupPD", number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("state") %>%
    rename(
      transformed_effect = logFC,
      average_transformed_proportion = AveExpr,
      p_value = P.Value,
      FDR = adj.P.Val
    ) %>%
    select(state, transformed_effect, average_transformed_proportion, t, p_value, FDR, B)
  means <- comp %>%
    group_by(state, group) %>%
    summarise(
      mean_proportion = mean(proportion),
      median_proportion = median(proportion),
      n_samples = n_distinct(sample),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = group,
      values_from = c(mean_proportion, median_proportion, n_samples)
    )
  result %>%
    left_join(means, by = "state") %>%
    mutate(
      subset = subset_name,
      formula = formula_string,
      n_samples_total = nrow(samples)
    )
}

message("Loading integrated parent object...")
parent <- readRDS(input_rds)
stopifnot(all(c("cellType", "sample", "group", "dataset") %in% colnames(parent[[]])))

message("Subsetting fibroblast lineage from the explicit-Harmony parent object...")
fibro <- subset(parent, subset = cellType == "Fibroblasts")
rm(parent)
invisible(gc())
DefaultAssay(fibro) <- "RNA"
fibro[["RNA"]] <- JoinLayers(fibro[["RNA"]])
fibro$cohort <- cohort_from_sample(as.character(fibro$sample))
fibro$cohort <- factor(fibro$cohort)
fibro$dataset <- factor(fibro$dataset)
fibro$group <- factor(fibro$group, levels = c("NC", "PD"))

if (file.exists(old_fibro_rds)) {
  message("Adding old fibroblast labels for sensitivity comparison...")
  old_fibro <- readRDS(old_fibro_rds)
  old_subtype <- as.character(old_fibro$subtype)
  names(old_subtype) <- colnames(old_fibro)
  fibro$old_subtype <- unname(old_subtype[colnames(fibro)])
  rm(old_fibro, old_subtype)
}

input_summary <- fibro[[]] %>%
  dplyr::count(dataset, cohort, group, sample, name = "n_fibro_candidates")
write.csv(input_summary, file.path(table_dir, "01_input_cells_by_sample.csv"), row.names = FALSE)
message("Fibroblast candidates: ", ncol(fibro), " cells from ", n_distinct(fibro$sample), " samples")
print(addmargins(table(fibro$dataset, fibro$group)))

message("Recomputing normalization, HVGs, cell-cycle/stress scores and PCA within fibroblasts...")
fibro <- NormalizeData(fibro, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
fibro <- FindVariableFeatures(fibro, selection.method = "vst", nfeatures = 2500, verbose = FALSE)

s_genes <- intersect(cc.genes.updated.2019$s.genes, rownames(fibro))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(fibro))
fibro <- CellCycleScoring(
  fibro,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)

programs <- list(
  inflammatory = c(
    "ICAM1", "CCL2", "CCL7", "CCL8", "CXCL1", "CXCL2", "CXCL3",
    "CXCL6", "CXCL8", "CXCL10", "CXCL12", "CXCL13", "IL6", "PTGS2", "NFKBIA"
  ),
  matrix_activated = c(
    "POSTN", "COL11A1", "CTHRC1", "TNC", "SERPINE1", "FAP", "SPARC",
    "COL1A1", "COL1A2", "COL3A1", "FN1", "THBS2", "ACTA2", "TAGLN", "MYL9"
  ),
  pi16_col15 = c(
    "PI16", "COL15A1", "MFAP5", "DPT", "CFD", "COL14A1", "C7", "C3"
  ),
  adventitial = c(
    "APOD", "PTGDS", "CFD", "C7", "ABCA8", "SFRP1", "COL14A1", "MGP"
  ),
  stress = c(
    "FOS", "JUN", "JUNB", "IER2", "IER3", "DNAJB1", "HSPA1A", "HSPA1B",
    "HSP90AA1", "ATF3", "DDIT3"
  ),
  myeloid_contamination = c("LST1", "TYROBP", "FCER1G", "PTPRC", "S100A8", "S100A9"),
  epithelial_contamination = c("EPCAM", "KRT8", "KRT14", "KRT18", "KRT19", "SFN"),
  endothelial_contamination = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "PLVAP"),
  fibroblast_identity = c("PDGFRA", "DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2")
)

program_gene_table <- bind_rows(lapply(names(programs), function(program_name) {
  data.frame(
    program = program_name,
    gene = programs[[program_name]],
    present = programs[[program_name]] %in% rownames(fibro)
  )
}))
write.csv(program_gene_table, file.path(table_dir, "02_independent_marker_programs.csv"), row.names = FALSE)

for (program_name in names(programs)) {
  present_genes <- intersect(programs[[program_name]], rownames(fibro))
  if (length(present_genes) < 3L) next
  raw_name <- paste0("program_", program_name)
  fibro <- AddModuleScore(
    fibro,
    features = list(present_genes),
    name = raw_name,
    seed = seed,
    search = FALSE
  )
  generated <- paste0(raw_name, "1")
  fibro[[raw_name]] <- fibro[[generated, drop = TRUE]]
}

fibro <- ScaleData(fibro, features = VariableFeatures(fibro), verbose = FALSE)
fibro <- RunPCA(
  fibro,
  features = VariableFeatures(fibro),
  npcs = 30,
  seed.use = seed,
  verbose = FALSE
)

message("Running fibroblast-specific Harmony by dataset...")
fibro <- RunHarmony(
  object = fibro,
  group.by.vars = "dataset",
  reduction = "pca",
  dims.use = 1:30,
  theta = 2,
  lambda = 1,
  sigma = 0.1,
  max_iter = 20,
  early_stop = TRUE,
  reduction.save = "harmony_fibro_dataset",
  plot_convergence = FALSE,
  verbose = TRUE
)

message("Running cohort-corrected Harmony as a sensitivity analysis...")
fibro <- RunHarmony(
  object = fibro,
  group.by.vars = "cohort",
  reduction = "pca",
  dims.use = 1:30,
  theta = 2,
  lambda = 1,
  sigma = 0.1,
  max_iter = 20,
  early_stop = TRUE,
  reduction.save = "harmony_fibro_cohort",
  plot_convergence = FALSE,
  verbose = TRUE
)

message("Building UMAPs and graphs...")
fibro <- RunUMAP(
  fibro,
  reduction = "pca",
  dims = 1:20,
  reduction.name = "umap_fibro_pca",
  seed.use = seed,
  verbose = FALSE
)
fibro <- RunUMAP(
  fibro,
  reduction = "harmony_fibro_dataset",
  dims = 1:20,
  reduction.name = "umap_fibro_harmony",
  seed.use = seed,
  verbose = FALSE
)
fibro <- FindNeighbors(
  fibro,
  reduction = "harmony_fibro_dataset",
  dims = 1:20,
  k.param = 20,
  graph.name = c("fibro_harmony_nn", "fibro_harmony_snn"),
  verbose = FALSE
)
fibro <- FindClusters(
  fibro,
  graph.name = "fibro_harmony_snn",
  resolution = c(0.2, 0.3, 0.4, 0.5),
  algorithm = 1,
  random.seed = seed,
  verbose = FALSE
)

primary_cluster_col <- "fibro_harmony_snn_res.0.3"
stopifnot(primary_cluster_col %in% colnames(fibro[[]]))
fibro$cluster_primary <- factor(fibro[[primary_cluster_col, drop = TRUE]])

fibro <- RunUMAP(
  fibro,
  reduction = "harmony_fibro_cohort",
  dims = 1:20,
  reduction.name = "umap_fibro_cohort",
  seed.use = seed,
  verbose = FALSE
)
fibro <- FindNeighbors(
  fibro,
  reduction = "harmony_fibro_cohort",
  dims = 1:20,
  k.param = 20,
  graph.name = c("fibro_cohort_nn", "fibro_cohort_snn"),
  verbose = FALSE
)
fibro <- FindClusters(
  fibro,
  graph.name = "fibro_cohort_snn",
  resolution = 0.3,
  algorithm = 1,
  random.seed = seed,
  cluster.name = "cluster_cohort_sensitivity",
  verbose = FALSE
)

resolution_cols <- grep("^fibro_harmony_snn_res\\.", colnames(fibro[[]]), value = TRUE)
resolution_ari <- expand_grid(a = resolution_cols, b = resolution_cols) %>%
  rowwise() %>%
  mutate(ARI = adjustedRandIndex(fibro[[a, drop = TRUE]], fibro[[b, drop = TRUE]])) %>%
  ungroup()
write.csv(resolution_ari, file.path(table_dir, "03_resolution_adjusted_rand_index.csv"), row.names = FALSE)

integration_ari <- data.frame(
  comparison = "dataset-Harmony cluster vs cohort-Harmony cluster",
  ARI = adjustedRandIndex(fibro$cluster_primary, fibro$cluster_cohort_sensitivity)
)
write.csv(integration_ari, file.path(table_dir, "04_integration_sensitivity_ARI.csv"), row.names = FALSE)

message("Ranking descriptive markers for new clusters...")
Idents(fibro) <- fibro$cluster_primary
markers <- FindAllMarkers(
  fibro,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25,
  max.cells.per.ident = 1000,
  random.seed = seed,
  verbose = FALSE
)
write.csv(markers, file.path(table_dir, "05_new_cluster_markers_all.csv"), row.names = FALSE)

specific_markers <- markers %>%
  mutate(pct_gap = pct.1 - pct.2) %>%
  filter(avg_log2FC >= 0.5, pct_gap >= 0.15) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30, with_ties = FALSE) %>%
  ungroup()
write.csv(specific_markers, file.path(table_dir, "06_new_cluster_markers_specific_top30.csv"), row.names = FALSE)

message("Summarizing cluster support and independent marker programs...")
metadata <- fibro[[]] %>% mutate(cell = rownames(fibro[[]]))
program_cols <- grep("^program_", colnames(metadata), value = TRUE)
for (program_col in program_cols) {
  metadata[[paste0(program_col, "_z")]] <- as.numeric(scale(metadata[[program_col]]))
}
program_z_cols <- paste0(program_cols, "_z")

cluster_support <- metadata %>%
  group_by(cluster_primary) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(sample),
    n_datasets = n_distinct(dataset),
    n_cohorts = n_distinct(cohort),
    max_dataset_fraction = max(prop.table(table(dataset))),
    max_cohort_fraction = max(prop.table(table(cohort))),
    max_sample_fraction = max(prop.table(table(sample))),
    NC_fraction = mean(group == "NC"),
    PD_fraction = mean(group == "PD"),
    median_percent_mt = median(percent.mt),
    median_nFeature_RNA = median(nFeature_RNA),
    median_S_score = median(S.Score),
    median_G2M_score = median(G2M.Score),
    .groups = "drop"
  )

cluster_program_means <- metadata %>%
  group_by(cluster_primary) %>%
  summarise(across(all_of(program_z_cols), mean), .groups = "drop")

biological_programs <- intersect(
  paste0("program_", c("inflammatory", "matrix_activated", "pi16_col15", "adventitial"), "_z"),
  names(cluster_program_means)
)
bio_matrix <- as.matrix(cluster_program_means[, biological_programs, drop = FALSE])
dominant_index <- max.col(bio_matrix, ties.method = "first")
dominant_program <- gsub("^program_|_z$", "", biological_programs[dominant_index])
second_score <- apply(bio_matrix, 1, function(x) sort(x, decreasing = TRUE)[min(2L, length(x))])
cluster_program_means$dominant_program <- dominant_program
cluster_program_means$program_margin <- apply(bio_matrix, 1, max) - second_score

label_dictionary <- c(
  inflammatory = "Inflammatory-like",
  matrix_activated = "Matrix-activated-like",
  pi16_col15 = "PI16/COL15A1-like",
  adventitial = "Adventitial-like"
)
cluster_program_means$preliminary_state <- unname(label_dictionary[cluster_program_means$dominant_program])

cluster_summary <- cluster_support %>%
  left_join(cluster_program_means, by = "cluster_primary") %>%
  mutate(
    annotation_confidence = case_when(
      max_sample_fraction >= 0.50 ~ "Low: sample-dominated",
      max_dataset_fraction >= 0.90 ~ "Low: dataset-limited",
      program_margin < 0.20 ~ "Low: mixed programs",
      TRUE ~ "Moderate"
    )
  )
write.csv(cluster_summary, file.path(table_dir, "07_new_cluster_support_and_programs.csv"), row.names = FALSE)

state_map <- setNames(cluster_summary$preliminary_state, as.character(cluster_summary$cluster_primary))
confidence_map <- setNames(cluster_summary$annotation_confidence, as.character(cluster_summary$cluster_primary))
fibro$fibro_state_preliminary <- factor(unname(state_map[as.character(fibro$cluster_primary)]))
fibro$annotation_confidence <- unname(confidence_map[as.character(fibro$cluster_primary)])
if (anyNA(fibro$fibro_state_preliminary)) {
  stop("Preliminary fibroblast-state mapping generated NA labels.")
}

metadata <- fibro[[]] %>% mutate(cell = rownames(fibro[[]]))
cluster_by_dataset <- metadata %>%
  dplyr::count(cluster_primary, dataset, name = "n") %>%
  group_by(cluster_primary) %>%
  mutate(fraction_within_cluster = n / sum(n)) %>%
  ungroup()
write.csv(cluster_by_dataset, file.path(table_dir, "08_cluster_by_dataset.csv"), row.names = FALSE)

cluster_by_sample <- metadata %>%
  dplyr::count(cluster_primary, dataset, cohort, group, sample, name = "n")
write.csv(cluster_by_sample, file.path(table_dir, "09_cluster_by_sample.csv"), row.names = FALSE)

message("Building sample-level composition tables and adjusted statistics...")
saveRDS(
  fibro,
  file.path(object_dir, "fibroblast_harmony_reanalysis_checkpoint.rds"),
  compress = FALSE
)
sample_composition <- composition_long(metadata, "fibro_state_preliminary", min_cells = 0L)
write.csv(sample_composition, file.path(table_dir, "10_sample_level_state_composition_all.csv"), row.names = FALSE)

stats_all <- bind_rows(lapply(c(50L, 100L), function(threshold) {
  comp <- sample_composition %>% filter(total_fibro >= threshold)
  run_propeller_style(
    comp,
    subset_name = paste0("All cohorts; minimum ", threshold, " fibroblasts/sample"),
    formula_string = "~ cohort + group"
  )
}))

stats_gse164241 <- bind_rows(lapply(c(50L, 100L), function(threshold) {
  comp <- sample_composition %>%
    filter(dataset == "GSE164241", total_fibro >= threshold)
  run_propeller_style(
    comp,
    subset_name = paste0("GSE164241; minimum ", threshold, " fibroblasts/sample"),
    formula_string = "~ cohort + group"
  )
}))

composition_statistics <- bind_rows(stats_all, stats_gse164241)
write.csv(
  composition_statistics,
  file.path(table_dir, "11_propeller_style_limma_composition_statistics.csv"),
  row.names = FALSE
)

message("Creating publication-oriented figures...")
dataset_colors <- c(GSE152042 = "#E76F51", GSE164241 = "#2A9D8F", GSE171213 = "#457B9D")
group_colors <- c(NC = "#5AB4AC", PD = "#D6604D")

p_before <- DimPlot(
  fibro,
  reduction = "umap_fibro_pca",
  group.by = "dataset",
  raster = TRUE,
  cols = dataset_colors
) +
  ggtitle("Before fibroblast-specific integration: PCA") +
  theme_classic(base_size = 11)

p_after <- DimPlot(
  fibro,
  reduction = "umap_fibro_harmony",
  group.by = "dataset",
  raster = TRUE,
  cols = dataset_colors
) +
  ggtitle("After fibroblast-specific Harmony: dataset") +
  theme_classic(base_size = 11)
save_plot("01_fibro_dataset_before_after.pdf", p_before + p_after, 13, 5.5)

p_dataset_split <- DimPlot(
  fibro,
  reduction = "umap_fibro_harmony",
  group.by = "cluster_primary",
  split.by = "dataset",
  raster = TRUE,
  ncol = 3
) +
  ggtitle("New clusters shown separately in each dataset") &
  theme_classic(base_size = 10)
save_plot("02_fibro_clusters_split_by_dataset.pdf", p_dataset_split, 16, 5.5)

p_cluster <- DimPlot(
  fibro,
  reduction = "umap_fibro_harmony",
  group.by = "cluster_primary",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("Fibroblast-specific Harmony clusters (resolution 0.3)") +
  theme_classic(base_size = 11)
p_group <- DimPlot(
  fibro,
  reduction = "umap_fibro_harmony",
  group.by = "group",
  cols = group_colors,
  raster = TRUE
) +
  ggtitle("Disease group") +
  theme_classic(base_size = 11)
save_plot("03_fibro_clusters_and_group.pdf", p_cluster + p_group, 13, 5.5)

p_old <- DimPlot(
  fibro,
  reduction = "umap_fibro_harmony",
  group.by = "old_subtype",
  na.value = "grey85",
  raster = TRUE
) +
  ggtitle("Previous fibroblast annotation") +
  theme_classic(base_size = 10)
p_new <- DimPlot(
  fibro,
  reduction = "umap_fibro_harmony",
  group.by = "fibro_state_preliminary",
  raster = TRUE
) +
  ggtitle("Independent-program preliminary annotation") +
  theme_classic(base_size = 10)
save_plot("04_previous_vs_preliminary_annotation.pdf", p_old + p_new, 14, 5.5)

marker_display <- unique(c(
  "PDGFRA", "DCN", "LUM",
  "ICAM1", "CCL2", "CXCL13", "CXCL12",
  "POSTN", "CTHRC1", "COL11A1", "TNC", "SERPINE1", "ACTA2",
  "PI16", "COL15A1", "MFAP5", "DPT",
  "APOD", "PTGDS", "CFD", "COL14A1",
  "FOS", "JUN", "IER3", "S100A8", "S100A9", "HBB", "EPCAM", "VWF"
))
marker_display <- intersect(marker_display, rownames(fibro))
p_dot <- DotPlot(fibro, features = marker_display, group.by = "cluster_primary") +
  RotatedAxis() +
  labs(
    title = "Independent fibroblast-state, stress and contamination markers",
    x = NULL,
    y = "New Harmony cluster"
  ) +
  theme_classic(base_size = 9)
save_plot("05_independent_marker_dotplot.pdf", p_dot, 15, 6.5)

sample_order <- sample_composition %>%
  distinct(sample, dataset, cohort, group) %>%
  arrange(dataset, cohort, group, sample) %>%
  pull(sample)
sample_plot_data <- sample_composition %>%
  mutate(sample = factor(sample, levels = sample_order))
p_sample <- ggplot(sample_plot_data, aes(x = sample, y = proportion, fill = state)) +
  geom_col(width = 0.9, color = "white", linewidth = 0.1) +
  facet_grid(~ dataset + group, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Fibroblast-state composition by biological sample",
    subtitle = "Each bar is one sample; samples, not cells, are the replicates",
    x = NULL,
    y = "Within-sample proportion",
    fill = "Preliminary state"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right")
save_plot("06_sample_level_state_composition.pdf", p_sample, 17, 6.5)

sample_valid <- sample_composition %>% filter(total_fibro >= 50)
p_group_state <- ggplot(
  sample_valid,
  aes(x = group, y = proportion, color = group)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.15) +
  geom_point(aes(shape = cohort), position = position_jitter(width = 0.12), size = 2, alpha = 0.85) +
  facet_wrap(~ state, scales = "free_y", ncol = 2) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Sample-level fibroblast-state proportions",
    subtitle = "Samples with at least 50 fibroblast candidates; shape denotes cohort",
    x = NULL,
    y = "Within-sample proportion",
    color = NULL,
    shape = "Cohort"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")
save_plot("07_state_proportions_by_group_and_cohort.pdf", p_group_state, 11, 8)

p_support <- ggplot(
  cluster_by_dataset,
  aes(x = dataset, y = cluster_primary, fill = fraction_within_cluster)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(fraction_within_cluster, accuracy = 0.1)), size = 3) +
  scale_fill_viridis_c(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Dataset support for each new fibroblast cluster",
    x = NULL,
    y = "New Harmony cluster",
    fill = "Fraction"
  ) +
  theme_classic(base_size = 11)
save_plot("08_cluster_dataset_support_heatmap.pdf", p_support, 8, 5.5)

message("Saving revised fibroblast object...")
saveRDS(
  fibro,
  file.path(object_dir, "fibroblast_harmony_reanalysis_preliminary.rds"),
  compress = FALSE
)

parameters <- data.frame(
  parameter = c(
    "seed", "input_parent", "input_fibroblast_definition", "n_fibroblast_candidates",
    "normalization", "HVGs", "PCA_dimensions", "Harmony_primary_variable",
    "Harmony_sensitivity_variable", "Harmony_theta", "Harmony_lambda", "Harmony_sigma",
    "Harmony_max_iter", "neighbor_dimensions", "k.param", "cluster_algorithm",
    "cluster_resolutions", "primary_resolution", "composition_unit", "composition_min_cells",
    "composition_model"
  ),
  value = c(
    seed, input_rds, "parent cellType == Fibroblasts", ncol(fibro),
    "LogNormalize scale.factor=10000", 2500, "1:30", "dataset",
    "cohort (sensitivity only)", 2, 1, 0.1, 20, "1:20", 20, "Louvain algorithm=1",
    "0.2,0.3,0.4,0.5", 0.3, "biological sample", "50 and 100",
    "arcsin-sqrt proportion + limma, design ~ cohort + group"
  )
)
write.csv(parameters, file.path(table_dir, "12_analysis_parameters.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
message("Primary clusters: ")
print(table(fibro$cluster_primary))
message("Preliminary states: ")
print(table(fibro$fibro_state_preliminary))
message("Integration sensitivity ARI: ", integration_ari$ARI)
message("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
