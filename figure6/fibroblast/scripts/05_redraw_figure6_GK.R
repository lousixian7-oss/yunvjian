args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1L) args[[1]] else "G:/1Yunvjian/0a26.7.7singlecell/8.10"
output_dir <- if (length(args) >= 2L) args[[2]] else file.path(root_dir, "fibroblast_revision_final")

object_file <- file.path(output_dir, "results/objects/fibroblast_revision_final.rds")
figure_dir <- file.path(output_dir, "results/figures/figure6_revised_GK")
table_dir <- file.path(output_dir, "results/tables/figure6_revised_GK")
log_dir <- file.path(output_dir, "logs")
invisible(lapply(c(figure_dir, table_dir, log_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

log_con <- file(file.path(log_dir, "05_redraw_figure6_GK.log"), "wt")
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
})

save_plot <- function(filename, plot, width, height) {
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  ggsave(pdf_path, plot, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(png_path, plot, width = width, height = height,
         dpi = 300, limitsize = FALSE)
}

message("Loading final fibroblast object...")
fibro <- readRDS(object_file)
DefaultAssay(fibro) <- "RNA"
fibro[["RNA"]] <- JoinLayers(fibro[["RNA"]])

defined_states <- c(
  "Inflammatory_Fib", "Activated_Fib", "ECM_Fib",
  "Adventitial_Fib", "PI16_Fib"
)
# Cluster 8 was deliberately left unresolved in the audit object.  For the final
# five-state display it is assigned to ECM_Fib because 93/132 cells have an
# ECM_Fib 30-nearest-neighbour majority in Harmony space and the ECM centroid is
# the closest of the five retained states.  This affects only 0.82% of cells.
fibro$state_final_audit <- as.character(fibro$state_final)
fibro$state_final <- as.character(fibro$state_final)
fibro$state_final[fibro$state_final == "Unassigned_Fib"] <- "ECM_Fib"
fibro$state_final <- factor(fibro$state_final, levels = defined_states)
all_states <- defined_states
write.csv(
  data.frame(
    original_state = "Unassigned_Fib", final_state = "ECM_Fib", n_cells = 132,
    fraction_all_fibroblasts = 132 / ncol(fibro),
    assignment_rule = paste(
      "Harmony 30-NN majority (93/132 ECM); closest state centroid = ECM;",
      "rare mixed cluster retained in audit metadata as state_final_audit"
    )
  ),
  file.path(table_dir, "G_rare_cluster_final_assignment.csv"), row.names = FALSE
)
saveRDS(fibro, file.path(output_dir,
                         "results/objects/fibroblast_revision_final_no_unassigned.rds"))

# Preserve the palette used in the original Figure 6 G-I panels.
fib_colors <- c(
  ECM_Fib = "#F8766D",
  PI16_Fib = "#CCCC00",
  Inflammatory_Fib = "#00BFC4",
  Activated_Fib = "#3399FF",
  Adventitial_Fib = "#C77CFF"
)
# Same hue family as the original Figure 6 palette, darkened slightly for the
# dense UMAP so individual cells remain visible after manuscript down-scaling.
fib_umap_colors <- c(
  ECM_Fib = "#E85C50",
  PI16_Fib = "#B3B300",
  Inflammatory_Fib = "#00A5AE",
  Activated_Fib = "#247FD1",
  Adventitial_Fib = "#AE63DB"
)
group_colors <- c(NC = "#7DD5D8", PD = "#F7B6AD")

base_theme <- theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.title = element_text(face = "plain", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "right"
  )

## G. Final fibroblast-state UMAP after cleaning and Harmony re-clustering.
p_g <- DimPlot(
  fibro,
  reduction = "umap_fibro_final",
  group.by = "state_final",
  cols = fib_umap_colors,
  label = TRUE,
  label.size = 4,
  repel = FALSE,
  pt.size = 0.90,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  labs(x = "umap_1", y = "umap_2", color = NULL, title = NULL) +
  ggtitle(NULL) +
  base_theme +
  theme(aspect.ratio = 1)

## H. Fixed marker annotation, matching the original grouped DotPlot layout.
marker_list <- list(
  Inflammatory_Fib = c("CCL2", "CXCL13"),
  Activated_Fib = c("SERPINE1", "TNC", "POSTN", "COL11A1"),
  ECM_Fib = c("COL1A1", "COL3A1", "LUM"),
  Adventitial_Fib = c("CFD", "PTGDS", "APOD"),
  PI16_Fib = c("PI16", "MFAP5", "COL15A1", "DPT")
)
marker_list <- lapply(marker_list, intersect, y = rownames(fibro))
fibro_defined <- subset(fibro, subset = state_final %in% defined_states)
fibro_defined$state_final <- droplevels(fibro_defined$state_final)
Idents(fibro_defined) <- fibro_defined$state_final
p_h <- DotPlot(
  fibro_defined,
  features = marker_list,
  group.by = "state_final",
  assay = "RNA",
  dot.scale = 8,
  cols = c("grey85", "blue")
) +
  RotatedAxis() +
  labs(x = "Features", y = "Identity", color = "Average Expression",
       size = "Percent Expressed") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 14),
    strip.text.x = element_text(size = 12),
    panel.spacing.x = unit(0.25, "lines"),
    legend.position = "right"
  )

## I. Replace pooled-cell composition with equal-weighted sample summaries.
sample_meta <- fibro[[]] %>%
  transmute(sample = as.character(sample), dataset = as.character(dataset),
            group = as.character(group)) %>%
  distinct()
composition <- fibro[[]] %>%
  transmute(sample = as.character(sample), state = as.character(state_final)) %>%
  count(sample, state, name = "n_state") %>%
  right_join(expand_grid(sample = sample_meta$sample, state = defined_states),
             by = c("sample", "state")) %>%
  mutate(n_state = replace_na(n_state, 0L)) %>%
  left_join(sample_meta, by = "sample") %>%
  group_by(sample) %>%
  mutate(total_fibro = sum(n_state), fraction = n_state / total_fibro) %>%
  ungroup() %>%
  mutate(
    eligible = total_fibro >= 50,
    state = factor(state, levels = defined_states),
    group = factor(group, levels = c("NC", "PD"))
  )

dataset_inclusion <- composition %>%
  distinct(sample, dataset, group, total_fibro, eligible) %>%
  arrange(dataset, group, sample)
write.csv(dataset_inclusion,
          file.path(table_dir, "I_all_GSM_dataset_inclusion.csv"), row.names = FALSE)
write.csv(composition,
          file.path(table_dir, "I_all_GSM_sample_state_composition.csv"), row.names = FALSE)

group_mean_composition <- composition %>%
  filter(eligible) %>%
  group_by(group, state) %>%
  summarise(
    mean_fraction = mean(fraction),
    n_samples = n_distinct(sample),
    .groups = "drop"
  )
write.csv(group_mean_composition,
          file.path(table_dir, "I_group_mean_sample_level_composition.csv"), row.names = FALSE)

p_i_stack <- ggplot(group_mean_composition,
                    aes(x = group, y = mean_fraction, fill = state)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = fib_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
  labs(x = NULL, y = "Mean fibroblast proportion", fill = "subtype",
       title = "Equal-weighted eligible GSMs") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

activated_fraction <- composition %>%
  filter(state == "Activated_Fib")
activated_eligible <- activated_fraction %>% filter(eligible)
activated_stats <- read.csv(
  file.path(output_dir, "results/tables/18_sample_level_composition_statistics.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
) %>%
  filter(state == "Activated_Fib", subset == "Integrated eligible samples", minimum_cells == 50)
act_subtitle <- sprintf("dataset-adjusted FDR = %.3g", activated_stats$FDR[1])
p_i_box <- ggplot(activated_fraction,
                  aes(x = group, y = fraction)) +
  geom_boxplot(data = activated_eligible, aes(fill = group), width = 0.58,
               outlier.shape = NA, linewidth = 0.55) +
  geom_point(aes(color = dataset, shape = eligible),
             position = position_jitter(width = 0.10), size = 2.8, stroke = 0.8) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = ">=50 cells (model)",
                                `FALSE` = "<50 cells (display only)")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Activated fibroblast fraction",
       title = "All three GSEs shown", subtitle = act_subtitle,
       fill = "group", color = "dataset", shape = "sample use") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

p_i <- p_i_stack + p_i_box + plot_layout(widths = c(0.92, 1.08))

## J/K. First summarize within each GSM x state, then average samples equally.
expr <- FetchData(
  fibro,
  vars = c("group", "state_final", "sample", "dataset", "C3", "AGT"),
  layer = "data"
) %>%
  mutate(state_final = as.character(state_final)) %>%
  filter(state_final %in% defined_states)

sample_expression <- expr %>%
  group_by(sample, dataset, group, state_final) %>%
  summarise(
    n_cells = n(),
    C3_avg = mean(C3), C3_pct = mean(C3 > 0) * 100,
    AGT_avg = mean(AGT), AGT_pct = mean(AGT > 0) * 100,
    .groups = "drop"
  )
write.csv(sample_expression,
          file.path(table_dir, "JK_sample_level_C3_AGT_by_state.csv"), row.names = FALSE)

# A minimum of 10 cells is required for a state-specific GSM estimate.
bubble_summary <- sample_expression %>%
  filter(n_cells >= 10) %>%
  group_by(group, state_final) %>%
  summarise(
    C3_avg = mean(C3_avg), C3_pct = mean(C3_pct),
    AGT_avg = mean(AGT_avg), AGT_pct = mean(AGT_pct),
    n_samples = n_distinct(sample),
    .groups = "drop"
  ) %>%
  complete(
    group = c("NC", "PD"), state_final = defined_states,
    fill = list(C3_avg = 0, C3_pct = 0, AGT_avg = 0, AGT_pct = 0, n_samples = 0)
  ) %>%
  mutate(
    group = factor(group, levels = c("NC", "PD")),
    state_final = factor(state_final, levels = defined_states)
  )
write.csv(bubble_summary,
          file.path(table_dir, "JK_equal_weighted_C3_AGT_bubble_summary.csv"), row.names = FALSE)

p_j <- ggplot(bubble_summary, aes(x = group, y = state_final)) +
  geom_point(aes(size = C3_pct, color = C3_avg)) +
  scale_size(range = c(1, 12), name = "C3+ cells (%)") +
  scale_color_gradient(low = "#00BFC4", high = "salmon",
                       name = "Average\nC3 expression") +
  labs(x = NULL, y = "Fibroblast subtype",
       title = "C3 expression in fibroblast subtypes",
       subtitle = "Equal-weighted GSM summaries; >=10 state cells/sample") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

p_k <- ggplot(bubble_summary, aes(x = group, y = state_final)) +
  geom_point(aes(size = AGT_pct, color = AGT_avg)) +
  scale_size(range = c(1, 12), name = "AGT+ cells (%)") +
  scale_color_gradient(low = "#00BFC4", high = "salmon",
                       name = "Average\nAGT expression") +
  labs(x = NULL, y = "Fibroblast subtype",
       title = "AGT expression in fibroblast subtypes",
       subtitle = "Equal-weighted GSM summaries; >=10 state cells/sample") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

save_plot("Figure6G_final_fibroblast_UMAP", p_g, 7.0, 5.8)
save_plot("Figure6H_final_fibroblast_marker_dotplot", p_h, 11.2, 4.7)
save_plot("Figure6I_sample_level_fibroblast_composition", p_i, 10.8, 4.8)
save_plot("Figure6J_C3_by_final_fibroblast_state", p_j, 6.7, 5.2)
save_plot("Figure6K_AGT_by_final_fibroblast_state", p_k, 6.7, 5.2)

# Combined G-K sheet for direct replacement in the original Figure 6 layout.
p_g_tag <- p_g + labs(tag = "G") +
  theme(plot.tag = element_text(face = "bold", size = 24))
p_h_tag <- p_h + labs(tag = "H") +
  theme(plot.tag = element_text(face = "bold", size = 24))
p_i_stack_tag <- p_i_stack + labs(tag = "I") +
  theme(plot.tag = element_text(face = "bold", size = 24))
p_i_tag <- p_i_stack_tag + p_i_box + plot_layout(widths = c(0.92, 1.08))
p_j_tag <- p_j + labs(tag = "J") +
  theme(plot.tag = element_text(face = "bold", size = 24))
p_k_tag <- p_k + labs(tag = "K") +
  theme(plot.tag = element_text(face = "bold", size = 24))

combined_gk <- (p_g_tag + p_i_tag + plot_layout(widths = c(0.82, 1.18))) /
  p_h_tag /
  (p_j_tag + p_k_tag) +
  plot_layout(heights = c(1.0, 0.75, 0.92))
save_plot("Figure6_GK_revised_combined", combined_gk, 16.0, 14.8)

writeLines(c(
  "Figure 6 G-K fibroblast re-analysis workflow (final display)",
  "1. Extracted major-cell-type fibroblast candidates from the integrated object.",
  "2. Removed lineage-contaminated clusters/cells using fixed pericyte, immune, epithelial, plasma-cell and myeloid marker panels.",
  "3. Re-normalized, selected variable genes, scaled, ran PCA, and integrated the fibroblast subset with Harmony using dataset only; disease group was not used for clustering or integration.",
  "4. Re-clustered at tested resolutions and annotated five fibroblast states with fixed literature-based marker programs.",
  "5. The rare audit cluster (132/16,032 cells) was assigned to ECM_Fib for the five-state final figure because 93/132 cells had an ECM 30-NN majority and the ECM centroid was closest in Harmony space.",
  "6. Figure 6I summarizes proportions at the GSM/sample level. All three GSEs are displayed; samples with fewer than 50 fibroblasts are shown as open symbols but excluded from the dataset-adjusted inferential model.",
  "7. C3/AGT subtype bubbles first summarize expression within each GSM x state (minimum 10 cells) and then weight GSMs equally; they are descriptive localization panels, not cell-level differential-expression tests."
), file.path(figure_dir, "Figure6_GK_analysis_workflow.txt"))

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_05_figure6_GK.txt"))
message("Completed revised Figure 6 G-K panels.")
