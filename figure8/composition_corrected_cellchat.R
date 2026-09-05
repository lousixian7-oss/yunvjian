## Composition-corrected sensitivity analysis for Figure 8
## This script never overwrites the original CellChat objects or figures.

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

INPUT_RDS <- "yjsl_fullannotation_fibro.rds"
OUT_DIR <- Sys.getenv("CELLCHAT_OUT_DIR", unset = "composition_corrected")
CONDITION_COL <- "group"
CELLTYPE_COL <- "full_annotation"
CONDITIONS <- c("NC", "PD")
N_REPEATS <- as.integer(Sys.getenv("CELLCHAT_N_REPEATS", unset = "20"))
NBOOT <- as.integer(Sys.getenv("CELLCHAT_NBOOT", unset = "100"))
MAX_CELLS_PER_TYPE <- 1000L
MIN_CELLS_PER_TYPE <- 20L
BASE_SEED <- 20260822L

cell_keep <- c(
  "T cells", "B cells", "Macrophages", "Endothelial cells",
  "Plasma cells", "Epithelial cells", "Pericytes", "Mast cells",
  "Neutrophils", "pDCs", "PI16_Fib", "Activated_Fib", "ECM_Fib",
  "Inflammatory_Fib", "Adventitial_Fib"
)

target_lr <- c(
  "COL1A1_CD44",
  "FN1_ITGA5_ITGB1",
  "FN1_ITGA8_ITGB1",
  "FN1_ITGAV_ITGB1"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(INPUT_RDS))
obj <- readRDS(INPUT_RDS)
stopifnot(all(c(CONDITION_COL, CELLTYPE_COL) %in% colnames(obj@meta.data)))

meta0 <- obj@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  filter(
    .data[[CONDITION_COL]] %in% CONDITIONS,
    .data[[CELLTYPE_COL]] %in% cell_keep
  ) %>%
  mutate(
    condition = .data[[CONDITION_COL]],
    cell_type = .data[[CELLTYPE_COL]]
  )

cell_counts <- meta0 %>%
  count(condition, cell_type, name = "n_cells") %>%
  complete(condition = CONDITIONS, cell_type = cell_keep, fill = list(n_cells = 0L))
write.csv(cell_counts, file.path(OUT_DIR, "01_cell_counts_before_matching.csv"), row.names = FALSE)

## Only retain cell types present with enough cells in both conditions.
matched_n <- cell_counts %>%
  group_by(cell_type) %>%
  summarise(n_match = min(n_cells), .groups = "drop") %>%
  mutate(n_match = pmin(n_match, MAX_CELLS_PER_TYPE)) %>%
  filter(n_match >= MIN_CELLS_PER_TYPE)

if (!all(c("PI16_Fib", "Activated_Fib") %in% matched_n$cell_type)) {
  stop("PI16_Fib and Activated_Fib must each have at least MIN_CELLS_PER_TYPE cells in both groups.")
}
write.csv(matched_n, file.path(OUT_DIR, "02_matched_cells_per_type.csv"), row.names = FALSE)

run_cellchat <- function(seu, condition_name) {
  seu <- NormalizeData(seu, verbose = FALSE)
  data_input <- GetAssayData(seu, assay = "RNA", layer = "data")
  meta <- seu@meta.data
  meta$CellChat_annotation <- meta[[CELLTYPE_COL]]

  cc <- createCellChat(
    object = data_input,
    meta = meta,
    group.by = "CellChat_annotation"
  )
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(
    cc,
    type = "triMean",
    raw.use = TRUE,
    population.size = FALSE,    # estimates expression-based, per-cell-type signaling
    nboot = NBOOT,
    seed.use = BASE_SEED
  )
  cc <- filterCommunication(cc, min.cells = 10) # identical rule in NC and PD
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

sample_one_repeat <- function(repeat_id) {
  message(sprintf("[%s] Starting matched repeat %d/%d", format(Sys.time()), repeat_id, N_REPEATS))
  set.seed(BASE_SEED + repeat_id)

  selected <- meta0 %>%
    inner_join(matched_n, by = "cell_type") %>%
    group_by(condition, cell_type) %>%
    group_modify(~ slice_sample(.x, n = unique(.x$n_match))) %>%
    ungroup()

  ## Verify strict matching before inference.
  check <- selected %>% count(condition, cell_type) %>%
    pivot_wider(names_from = condition, values_from = n)
  stopifnot(all(check$NC == check$PD))

  out <- vector("list", length(CONDITIONS))
  names(out) <- CONDITIONS
  for (cond in CONDITIONS) {
    cells <- selected$barcode[selected$condition == cond]
    seu <- CreateSeuratObject(
      counts = GetAssayData(obj, assay = "RNA", layer = "counts")[, cells, drop = FALSE],
      meta.data = obj@meta.data[cells, , drop = FALSE]
    )
    out[[cond]] <- run_cellchat(seu, cond)
  }

  global <- bind_rows(lapply(CONDITIONS, function(cond) {
    cc <- out[[cond]]
    data.frame(
      repeat_id = repeat_id,
      condition = cond,
      total_interactions = sum(cc@net$count),
      total_weight = sum(cc@net$weight)
    )
  }))

  lr <- bind_rows(lapply(CONDITIONS, function(cond) {
    subsetCommunication(
      out[[cond]],
      sources.use = "PI16_Fib",
      targets.use = "Activated_Fib"
    ) %>%
      transmute(
        repeat_id = repeat_id,
        condition = cond,
        interaction_name,
        ligand,
        receptor,
        pathway_name,
        prob,
        within_condition_p = pval
      )
  })) %>%
    filter(interaction_name %in% target_lr)

  ## Save one representative pair of objects for revised panels; do not save all repeats.
  if (repeat_id == 1L) {
    saveRDS(out$NC, file.path(OUT_DIR, "CellChat_NC_matched_populationFALSE.rds"))
    saveRDS(out$PD, file.path(OUT_DIR, "CellChat_PD_matched_populationFALSE.rds"))
  }

  message(sprintf("[%s] Finished matched repeat %d/%d", format(Sys.time()), repeat_id, N_REPEATS))
  ans <- list(global = global, lr = lr)
  rm(out)
  invisible(gc())
  ans
}

res <- map(seq_len(N_REPEATS), sample_one_repeat)
global_all <- bind_rows(map(res, "global"))
lr_all <- bind_rows(map(res, "lr"))

write.csv(global_all, file.path(OUT_DIR, "03_global_metrics_all_resamples.csv"), row.names = FALSE)
write.csv(lr_all, file.path(OUT_DIR, "04_target_LR_all_resamples.csv"), row.names = FALSE)

global_delta <- global_all %>%
  pivot_wider(names_from = condition, values_from = c(total_interactions, total_weight)) %>%
  mutate(
    delta_count = total_interactions_PD - total_interactions_NC,
    delta_weight = total_weight_PD - total_weight_NC
  )
write.csv(global_delta, file.path(OUT_DIR, "05_global_PD_minus_NC.csv"), row.names = FALSE)

lr_delta <- lr_all %>%
  select(repeat_id, condition, interaction_name, ligand, receptor, pathway_name, prob) %>%
  complete(
    repeat_id,
    condition = CONDITIONS,
    nesting(interaction_name, ligand, receptor, pathway_name),
    fill = list(prob = 0)
  ) %>%
  pivot_wider(names_from = condition, values_from = prob, values_fill = 0) %>%
  mutate(delta_prob = PD - NC) %>%
  group_by(interaction_name, ligand, receptor, pathway_name) %>%
  summarise(
    median_delta = median(delta_prob),
    q025 = quantile(delta_prob, 0.025),
    q975 = quantile(delta_prob, 0.975),
    fraction_PD_higher = mean(delta_prob > 0),
    n_resamples = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(fraction_PD_higher), desc(median_delta))
write.csv(lr_delta, file.path(OUT_DIR, "06_target_LR_stability_summary.csv"), row.names = FALSE)

p_global <- global_delta %>%
  select(repeat_id, delta_count, delta_weight) %>%
  pivot_longer(-repeat_id, names_to = "metric", values_to = "PD_minus_NC") %>%
  ggplot(aes(x = metric, y = PD_minus_NC, fill = metric)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.5) +
  scale_x_discrete(labels = c(delta_count = "Interaction count", delta_weight = "Interaction weight")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "PD - NC (composition-matched resamples)")
ggsave(file.path(OUT_DIR, "Figure8A_composition_sensitivity.pdf"), p_global, width = 5, height = 4)

p_lr <- lr_delta %>%
  ggplot(aes(x = reorder(interaction_name, median_delta), y = median_delta)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.15) +
  geom_point(aes(colour = fraction_PD_higher), size = 3) +
  coord_flip() +
  scale_colour_viridis_c(limits = c(0, 1), name = "Fraction PD > NC") +
  theme_classic() +
  labs(x = "PI16_Fib to Activated_Fib L-R pair", y = "Median communication-probability difference")
ggsave(file.path(OUT_DIR, "Figure8DE_target_LR_stability.pdf"), p_lr, width = 7, height = 4.5)

message("Finished. Outputs are in: ", normalizePath(OUT_DIR))
message("Important: resampling intervals measure cell-sampling stability, not biological-replicate significance.")
