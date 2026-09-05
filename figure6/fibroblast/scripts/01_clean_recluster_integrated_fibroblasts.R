args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1L) args[[1]] else "G:/1Yunvjian/0a26.7.7singlecell/8.10"
output_dir <- if (length(args) >= 2L) args[[2]] else file.path(root_dir, "fibroblast_revision_final")
input_rds <- file.path(output_dir, "results/objects/fibroblast_preliminary_qc_checkpoint.rds")

figure_dir <- file.path(output_dir, "results/figures")
supp_dir <- file.path(figure_dir, "supplement")
table_dir <- file.path(output_dir, "results/tables")
object_dir <- file.path(output_dir, "results/objects")
log_dir <- file.path(output_dir, "logs")
invisible(lapply(c(figure_dir, supp_dir, table_dir, object_dir, log_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

log_con <- file(file.path(log_dir, "01_clean_recluster_integrated_fibroblasts.log"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_con)
}, add = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(mclust)
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
  generated <- paste0(score_name, "1")
  object[[score_name]] <- object[[generated, drop = TRUE]]
  object
}

cramers_v <- function(tab) {
  test <- suppressWarnings(chisq.test(tab, correct = FALSE))
  sqrt(as.numeric(test$statistic) /
         (sum(tab) * min(nrow(tab) - 1L, ncol(tab) - 1L)))
}

message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("Loading: ", input_rds)
fibro_before <- readRDS(input_rds)
DefaultAssay(fibro_before) <- "RNA"
fibro_before[["RNA"]] <- JoinLayers(fibro_before[["RNA"]])

contamination_sets <- list(
  Pericyte_smooth_muscle = c("GUCY1A1", "GUCY1B1", "RGS5", "MCAM"),
  T_NK = c("CD3D", "TRAC", "NKG7", "GZMK"),
  Epithelial = c("KRT5", "KRT6A", "KRT13", "DSP"),
  Plasma = c("MZB1", "CD79A", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGKC"),
  Myeloid_ambiguous = c("CD68", "ACP5")
)

for (nm in names(contamination_sets)) {
  fibro_before <- add_score(fibro_before, contamination_sets[[nm]], paste0("contam_", nm))
}
fibro_before <- add_score(
  fibro_before,
  unique(unlist(contamination_sets, use.names = FALSE)),
  "contamination_score"
)

decision <- data.frame(
  preliminary_cluster = as.character(0:11),
  decision = c(
    "Keep", "Keep", "Keep", "Keep",
    "Exclude", "Keep", "Keep", "Keep",
    "Exclude", "Exclude", "Exclude", "Exclude"
  ),
  assigned_identity = c(
    rep("Fibroblast candidate", 4),
    "Pericyte/smooth muscle",
    rep("Fibroblast candidate", 3),
    "Myeloid/ambiguous sample-dominated",
    "T/NK", "Epithelial", "Plasma"
  ),
  evidence = c(
    "Fibroblast markers retained", "Fibroblast markers retained",
    "Fibroblast markers retained", "Fibroblast markers retained",
    "GUCY1A1/GUCY1B1",
    "Fibroblast inflammatory markers", "Fibroblast markers; sex-associated split",
    "PI16/MFAP5/DPT",
    "ACP5/CD68; weak fibroblast identity; one-sample dominance",
    "CD3D/TRAC/NKG7/GZMK", "KRT5/KRT6A/KRT13/DSP",
    "MZB1/CD79A/immunoglobulin genes"
  ),
  stringsAsFactors = FALSE
)
write.csv(decision, file.path(table_dir, "01_contamination_cluster_decisions.csv"), row.names = FALSE)

decision_map <- setNames(decision$assigned_identity, decision$preliminary_cluster)
keep_clusters <- decision$preliminary_cluster[decision$decision == "Keep"]
fibro_before$qc_lineage_call <- unname(decision_map[as.character(fibro_before$cluster_primary)])
fibro_before$qc_keep <- as.character(fibro_before$cluster_primary) %in% keep_clusters

contam_markers <- intersect(unique(unlist(contamination_sets, use.names = FALSE)), rownames(fibro_before))
p_contam_dot <- DotPlot(fibro_before, features = contam_markers, group.by = "cluster_primary") +
  RotatedAxis() +
  labs(title = "Lineage-contamination markers before fibroblast cleaning",
       x = NULL, y = "Preliminary cluster") +
  theme_classic(base_size = 9)
save_plot(file.path(supp_dir, "S01_contamination_marker_dotplot.pdf"), p_contam_dot, 12, 6)

p_before <- DimPlot(fibro_before, reduction = "umap_fibro_harmony",
                    group.by = "qc_lineage_call", raster = TRUE) +
  ggtitle("Before lineage-contamination removal") +
  theme_classic(base_size = 10)

clean_cells <- colnames(fibro_before)[fibro_before$qc_keep]
fibro <- subset(fibro_before, cells = clean_cells)
rm(clean_cells)

cleaning_summary <- data.frame(
  metric = c("fibroblast_candidates_before", "cells_retained", "cells_excluded",
             "retention_fraction"),
  value = c(ncol(fibro_before), ncol(fibro), ncol(fibro_before) - ncol(fibro),
            ncol(fibro) / ncol(fibro_before))
)
write.csv(cleaning_summary, file.path(table_dir, "02_contamination_cleaning_summary.csv"), row.names = FALSE)

excluded_counts <- fibro_before[[]] %>%
  filter(!qc_keep) %>%
  dplyr::count(cluster_primary, qc_lineage_call, dataset, group, name = "n_cells")
write.csv(excluded_counts, file.path(table_dir, "03_excluded_cells_by_cluster_dataset_group.csv"), row.names = FALSE)

message("Retained ", ncol(fibro), " / ", ncol(fibro_before), " candidate cells")
DefaultAssay(fibro) <- "RNA"
fibro[["RNA"]] <- JoinLayers(fibro[["RNA"]])
fibro$dataset <- factor(fibro$dataset)
fibro$group <- factor(fibro$group, levels = c("NC", "PD"))
fibro$cohort <- factor(fibro$cohort)

message("Recomputing normalized expression, HVGs and PCA after cleaning...")
fibro <- NormalizeData(fibro, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
fibro <- FindVariableFeatures(fibro, selection.method = "vst", nfeatures = 3500,
                              verbose = FALSE)

stress_genes <- c("FOS", "JUN", "IER3", "NR4A3")
technical_regex <- paste0(
  "^(MT-|RPL|RPS|MTRNR|HBA|HBB|IGH|IGK|IGL|TRAV|TRBV|TRGV|TRDV|",
  "MIR|SNOR|RP11-|AC[0-9]|AL[0-9])"
)
technical_hvgs <- unique(c(
  grep(technical_regex, VariableFeatures(fibro), value = TRUE),
  "XIST", "TSIX", "DDX3Y", "RPS4Y1", "KDM5D", "UTY", "EIF1AY",
  cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes,
  stress_genes,
  unique(unlist(contamination_sets, use.names = FALSE))
))
technical_hvgs <- intersect(technical_hvgs, VariableFeatures(fibro))
biological_hvgs <- setdiff(VariableFeatures(fibro), technical_hvgs)
VariableFeatures(fibro) <- biological_hvgs[seq_len(min(2500L, length(biological_hvgs)))]
write.csv(data.frame(gene = technical_hvgs),
          file.path(table_dir, "04_technical_HVGs_excluded_from_clustering.csv"),
          row.names = FALSE)

fibro <- ScaleData(fibro, features = VariableFeatures(fibro), verbose = FALSE)
fibro <- RunPCA(fibro, features = VariableFeatures(fibro), npcs = 30,
                seed.use = seed, verbose = FALSE)

message("Running Harmony using dataset only; disease group is not used...")
fibro <- RunHarmony(
  fibro,
  group.by.vars = "dataset",
  reduction = "pca",
  dims.use = 1:30,
  theta = 2,
  lambda = 1,
  sigma = 0.1,
  max_iter = 20,
  early_stop = TRUE,
  reduction.save = "harmony_fibro_final",
  plot_convergence = FALSE,
  verbose = TRUE
)

fibro <- RunUMAP(fibro, reduction = "pca", dims = 1:20,
                 reduction.name = "umap_fibro_before_final",
                 seed.use = seed, verbose = FALSE)
fibro <- RunUMAP(fibro, reduction = "harmony_fibro_final", dims = 1:20,
                 reduction.name = "umap_fibro_final",
                 seed.use = seed, verbose = FALSE)
fibro <- FindNeighbors(
  fibro, reduction = "harmony_fibro_final", dims = 1:20, k.param = 20,
  graph.name = c("fibro_final_nn", "fibro_final_snn"), verbose = FALSE
)
fibro <- FindClusters(
  fibro, graph.name = "fibro_final_snn",
  resolution = c(0.2, 0.3, 0.4, 0.5), algorithm = 1,
  random.seed = seed, verbose = FALSE
)
primary_col <- "fibro_final_snn_res.0.3"
stopifnot(primary_col %in% colnames(fibro[[]]))
fibro$cluster_final <- factor(fibro[[primary_col, drop = TRUE]])

p_after <- DimPlot(fibro, reduction = "umap_fibro_final",
                   group.by = "cluster_final", label = TRUE, repel = TRUE,
                   raster = TRUE) +
  ggtitle("After contamination removal and fibroblast-specific Harmony") +
  theme_classic(base_size = 10)
save_plot(file.path(supp_dir, "S02_fibroblast_cleaning_before_after_umap.pdf"),
          p_before + p_after, 14, 5.5)

resolution_cols <- grep("^fibro_final_snn_res\\.", colnames(fibro[[]]), value = TRUE)
resolution_ari <- expand_grid(a = resolution_cols, b = resolution_cols) %>%
  rowwise() %>%
  mutate(ARI = adjustedRandIndex(fibro[[a, drop = TRUE]], fibro[[b, drop = TRUE]])) %>%
  ungroup()
write.csv(resolution_ari, file.path(table_dir, "05_resolution_stability_ARI.csv"), row.names = FALSE)

message("Computing fixed fibroblast-state programs...")
state_sets <- list(
  PI16_fibroblast = c("PI16", "COL15A1", "MFAP5", "DPT"),
  Activated_fibroblast = c("POSTN", "TNC", "COL1A1", "COL11A1", "SERPINE1", "FAP"),
  ECM_fibroblast = c("COL1A1", "COL3A1", "COL5A1"),
  Adventitial_fibroblast = c("COL14A1", "APOD", "PTGDS", "CFD"),
  Myofibroblast = c("ACTA2", "TAGLN", "MYL9")
)
for (nm in names(state_sets)) fibro <- add_score(fibro, state_sets[[nm]], paste0("score_", nm))
fibro <- add_score(fibro, stress_genes, "score_Stress")
fibro <- add_score(fibro, unique(unlist(contamination_sets, use.names = FALSE)),
                   "score_Contamination")

marker_program_table <- bind_rows(lapply(names(state_sets), function(nm) {
  data.frame(program = nm, gene = state_sets[[nm]], present = state_sets[[nm]] %in% rownames(fibro))
}))
marker_program_table <- bind_rows(
  marker_program_table,
  data.frame(program = "Stress", gene = stress_genes, present = stress_genes %in% rownames(fibro)),
  data.frame(program = "Contamination", gene = unique(unlist(contamination_sets, use.names = FALSE)),
             present = unique(unlist(contamination_sets, use.names = FALSE)) %in% rownames(fibro))
)
write.csv(marker_program_table, file.path(table_dir, "06_fixed_marker_programs.csv"), row.names = FALSE)

Idents(fibro) <- fibro$cluster_final
markers <- FindAllMarkers(
  fibro, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
  min.pct = 0.10, logfc.threshold = 0.25, max.cells.per.ident = 1000,
  random.seed = seed, verbose = FALSE
) %>% mutate(pct_gap = pct.1 - pct.2)
write.csv(markers, file.path(table_dir, "07_final_cluster_markers_all.csv"), row.names = FALSE)
top_markers <- markers %>%
  filter(avg_log2FC >= 0.5, pct_gap >= 0.15) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 30, with_ties = FALSE) %>%
  ungroup()
write.csv(top_markers, file.path(table_dir, "08_final_cluster_markers_specific_top30.csv"), row.names = FALSE)

metadata <- fibro[[]]
state_score_cols <- paste0("score_", names(state_sets))
for (score_col in c(state_score_cols, "score_Stress", "score_Contamination")) {
  metadata[[paste0(score_col, "_z")]] <- as.numeric(scale(metadata[[score_col]]))
}

cluster_support <- metadata %>%
  group_by(cluster_final) %>%
  summarise(
    n_cells = n(), n_samples = n_distinct(sample), n_datasets = n_distinct(dataset),
    n_cohorts = n_distinct(cohort),
    max_dataset_fraction = max(prop.table(table(dataset))),
    max_sample_fraction = max(prop.table(table(sample))),
    NC_fraction = mean(group == "NC"), PD_fraction = mean(group == "PD"),
    median_percent_mt = median(percent.mt),
    median_nFeature_RNA = median(nFeature_RNA),
    across(ends_with("_z"), mean),
    .groups = "drop"
  )

bio_z_cols <- paste0(state_score_cols, "_z")
bio_matrix <- as.matrix(cluster_support[, bio_z_cols, drop = FALSE])
dominant_index <- max.col(bio_matrix, ties.method = "first")
dominant_state <- names(state_sets)[dominant_index]
second_score <- apply(bio_matrix, 1, function(x) sort(x, decreasing = TRUE)[min(2L, length(x))])
cluster_support$dominant_fixed_program <- dominant_state
cluster_support$program_margin <- apply(bio_matrix, 1, max) - second_score
cluster_support$annotation_flag <- case_when(
  cluster_support$score_Contamination_z > 0.75 ~ "Contamination-high",
  cluster_support$score_Stress_z > 0.75 ~ "Stress-high",
  cluster_support$program_margin < 0.15 ~ "Mixed-program",
  cluster_support$max_dataset_fraction > 0.95 ~ "Dataset-limited",
  TRUE ~ "Supported"
)
write.csv(cluster_support, file.path(table_dir, "09_cluster_fixed_program_summary.csv"), row.names = FALSE)

state_map <- setNames(cluster_support$dominant_fixed_program,
                      as.character(cluster_support$cluster_final))
flag_map <- setNames(cluster_support$annotation_flag,
                     as.character(cluster_support$cluster_final))
fibro$state_fixed_preliminary <- factor(unname(state_map[as.character(fibro$cluster_final)]),
                                        levels = names(state_sets))
fibro$state_annotation_flag <- unname(flag_map[as.character(fibro$cluster_final)])

dataset_state_table <- table(fibro$dataset, fibro$state_fixed_preliminary)
write.csv(as.data.frame.matrix(dataset_state_table),
          file.path(table_dir, "10_dataset_by_preliminary_state_contingency.csv"))
dataset_association <- data.frame(
  association = c("dataset x cluster_final", "dataset x preliminary_state"),
  cramers_v = c(
    cramers_v(table(fibro$dataset, fibro$cluster_final)),
    cramers_v(dataset_state_table)
  )
)
write.csv(dataset_association, file.path(table_dir, "11_dataset_association_cramers_v.csv"),
          row.names = FALSE)

heat_data <- as.data.frame(dataset_state_table) %>%
  rename(dataset = Var1, state = Var2, n = Freq) %>%
  group_by(state) %>%
  mutate(fraction_within_state = n / sum(n)) %>%
  ungroup()
p_heat <- ggplot(heat_data, aes(x = dataset, y = state, fill = fraction_within_state)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(fraction_within_state, accuracy = 0.1)), size = 3) +
  scale_fill_viridis_c(limits = c(0, 1), labels = scales::percent) +
  labs(title = "Dataset association of fixed-program fibroblast states",
       x = NULL, y = NULL, fill = "Within-state\nfraction") +
  theme_classic(base_size = 10)
save_plot(file.path(supp_dir, "S03_dataset_state_association_heatmap.pdf"), p_heat, 8, 5)

marker_display <- intersect(unique(unlist(state_sets, use.names = FALSE)), rownames(fibro))
p_marker <- DotPlot(fibro, features = marker_display, group.by = "cluster_final") +
  RotatedAxis() +
  labs(title = "Fixed fibroblast-state markers", x = NULL, y = "Final cluster") +
  theme_classic(base_size = 9)
save_plot(file.path(figure_dir, "F6B_fixed_marker_dotplot.pdf"), p_marker, 13, 6)

p_umap <- DimPlot(fibro, reduction = "umap_fibro_final",
                  group.by = "state_fixed_preliminary", raster = TRUE) +
  ggtitle("Integrated fibroblast states defined by fixed marker programs") +
  theme_classic(base_size = 10)
save_plot(file.path(figure_dir, "F6A_integrated_fibroblast_umap_preliminary.pdf"), p_umap, 7, 5.5)

feature_scores <- c(state_score_cols, "score_Stress", "score_Contamination")
p_features <- FeaturePlot(fibro, features = feature_scores,
                          reduction = "umap_fibro_final", ncol = 3,
                          order = TRUE, raster = TRUE) & theme_classic(base_size = 8)
save_plot(file.path(figure_dir, "F6_module_score_featureplots.pdf"), p_features, 13, 10)

score_long <- fibro[[]] %>%
  select(cluster_final, all_of(feature_scores)) %>%
  pivot_longer(cols = all_of(feature_scores), names_to = "program", values_to = "score")
p_violin <- ggplot(score_long, aes(x = cluster_final, y = score, fill = cluster_final)) +
  geom_violin(scale = "width", linewidth = 0.15) +
  facet_wrap(~ program, scales = "free_y", ncol = 3) +
  labs(title = "Fixed module scores across final clusters", x = "Final cluster", y = "Module score") +
  theme_classic(base_size = 9) + theme(legend.position = "none")
save_plot(file.path(figure_dir, "F6_module_score_violins.pdf"), p_violin, 13, 9)

message("Saving cleaned and reclustered fibroblast object...")
saveRDS(fibro, file.path(object_dir, "fibroblast_clean_harmony_fixed_programs_preliminary.rds"),
        compress = FALSE)

parameters <- data.frame(
  parameter = c(
    "seed", "input", "preliminary_resolution", "retained_clusters", "excluded_clusters",
    "normalization", "initial_HVGs", "final_HVG_max", "technical_HVG_filter",
    "Harmony_variable", "Harmony_theta", "Harmony_lambda", "Harmony_sigma",
    "Harmony_max_iter", "PCA_dims", "neighbor_dims", "cluster_resolutions",
    "primary_resolution", "disease_group_used_in_clustering"
  ),
  value = c(
    seed, input_rds, 0.3, paste(keep_clusters, collapse = ","),
    paste(decision$preliminary_cluster[decision$decision == "Exclude"], collapse = ","),
    "LogNormalize scale.factor=10000", 3500, 2500,
    "mitochondrial/ribosomal/IG/TCR/sex/cell-cycle/stress/contamination genes",
    "dataset", 2, 1, 0.1, 20, "1:30", "1:20", "0.2,0.3,0.4,0.5", 0.3, "No"
  )
)
write.csv(parameters, file.path(table_dir, "12_clean_recluster_parameters.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_clean_recluster.txt"))
message("Final preliminary cluster sizes:")
print(table(fibro$cluster_final))
message("Fixed-program preliminary states:")
print(table(fibro$state_fixed_preliminary))
message("Dataset association:")
print(dataset_association)
message("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
