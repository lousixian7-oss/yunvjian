args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1L) args[[1]] else "G:/1Yunvjian/0a26.7.7singlecell/8.10"
output_dir <- if (length(args) >= 2L) args[[2]] else file.path(root_dir, "fibroblast_revision_final")

object_file <- file.path(output_dir, "results/objects/fibroblast_revision_final.rds")
figure_dir <- file.path(output_dir, "results/figures")
main_dir <- file.path(figure_dir, "main_figure6")
supp_dir <- file.path(figure_dir, "supplement")
table_dir <- file.path(output_dir, "results/tables")
log_dir <- file.path(output_dir, "logs")
invisible(lapply(c(main_dir, supp_dir, log_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

log_con <- file(file.path(log_dir, "03_assemble_final_figures.log"), "wt")
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

save_plot <- function(path, plot, width, height) {
  ggsave(path, plot, width = width, height = height, device = cairo_pdf, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", path), plot, width = width, height = height,
         dpi = 300, limitsize = FALSE)
}

message("Loading final fibroblast object...")
fibro <- readRDS(object_file)
DefaultAssay(fibro) <- "RNA"
fibro[["RNA"]] <- JoinLayers(fibro[["RNA"]])

state_levels <- c("PI16_Fib", "Activated_Fib", "ECM_Fib", "Adventitial_Fib",
                  "Inflammatory_Fib", "Unassigned_Fib")
state_labels <- c(
  PI16_Fib = "PI16 fibroblasts",
  Activated_Fib = "Activated fibroblasts",
  ECM_Fib = "ECM fibroblasts",
  Adventitial_Fib = "Adventitial fibroblasts",
  Inflammatory_Fib = "Inflammatory fibroblasts",
  Unassigned_Fib = "Unassigned fibroblasts"
)
state_colors <- c(
  PI16_Fib = "#D4A017", Activated_Fib = "#D55E00", ECM_Fib = "#7B61A8",
  Adventitial_Fib = "#009E73", Inflammatory_Fib = "#56B4E9", Unassigned_Fib = "#999999"
)
group_colors <- c(NC = "#4DBBD5", PD = "#E64B35")

fibro$state_final <- factor(as.character(fibro$state_final), levels = state_levels)
Idents(fibro) <- fibro$state_final

panel_theme <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8))

# A: final-state UMAP. Disease group was not used in clustering or Harmony.
p_a <- DimPlot(
  fibro, reduction = "umap_fibro_final", group.by = "state_final",
  cols = state_colors, label = TRUE, repel = TRUE, raster = TRUE
) +
  scale_color_manual(values = state_colors, labels = state_labels, drop = FALSE) +
  labs(title = "Integrated fibroblast states", subtitle = "Harmony correction used dataset only; disease group was not used",
       color = "Fibroblast state", x = "UMAP 1", y = "UMAP 2") +
  panel_theme

# B: fixed, pre-specified marker programs summarized by final state.
marker_genes <- c(
  "PI16", "COL15A1", "MFAP5", "DPT",
  "POSTN", "TNC", "COL11A1", "SERPINE1", "FAP",
  "COL1A1", "COL3A1", "COL5A1",
  "COL14A1", "APOD", "PTGDS", "CFD",
  "ACTA2", "TAGLN", "MYL9",
  "CXCL2", "CXCL6", "CXCL13", "IL6"
)
marker_genes <- intersect(marker_genes, rownames(fibro))
p_b <- DotPlot(fibro, features = marker_genes, group.by = "state_final",
               assay = "RNA", dot.scale = 5) +
  scale_y_discrete(labels = state_labels) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B") +
  RotatedAxis() +
  labs(title = "Fixed-marker annotation", x = NULL, y = NULL,
       color = "Scaled mean", size = "% expressed") +
  panel_theme + theme(axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1, size = 7.0),
                      axis.text.y = element_text(size = 8))

# C: sample-level composition; every bar is one GSM.
composition <- read.csv(file.path(table_dir, "16_sample_level_fibroblast_composition.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
composition$state <- factor(composition$state, levels = state_levels)
sample_order <- composition |>
  distinct(sample, dataset, cohort, group) |>
  arrange(dataset, cohort, group, sample) |>
  pull(sample)
composition$sample <- factor(composition$sample, levels = sample_order)
sample_totals <- composition |>
  distinct(sample, total_fibro) |>
  mutate(sample_label = paste0(sample, "\n(n=", total_fibro, ")"))
sample_labels <- setNames(sample_totals$sample_label, sample_totals$sample)
p_c <- ggplot(composition, aes(sample, fraction, fill = state)) +
  geom_col(width = 0.88, color = "white", linewidth = 0.1) +
  facet_grid(~ dataset + group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = state_colors, labels = state_labels, drop = FALSE) +
  scale_x_discrete(labels = sample_labels) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Fibroblast-state composition by biological sample",
       subtitle = "One bar per GSM; denominator = all qualified fibroblasts in that sample",
       x = NULL, y = "Within-sample fraction", fill = "Fibroblast state") +
  panel_theme +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 6.5),
        strip.text = element_text(face = "bold", size = 8), legend.position = "right")

# D: activated program is tested at the biological-sample level.
sample_program <- read.csv(file.path(table_dir, "17_sample_level_module_scores.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(total_fibro >= 50)
act_stats <- read.csv(file.path(table_dir, "19_sample_level_module_score_statistics.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE) |>
  filter(program == "score_Activated_fibroblast", subset == "Integrated eligible samples")
stat_subtitle <- sprintf("Sample-level means (n=%d); dataset-adjusted FDR=%.4f",
                         act_stats$n_samples_total[1], act_stats$FDR[1])
p_d <- ggplot(sample_program, aes(group, score_Activated_fibroblast, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.12, linewidth = 0.6) +
  geom_point(aes(shape = dataset), position = position_jitter(width = 0.10), size = 2.2, alpha = 0.9) +
  scale_color_manual(values = group_colors) +
  labs(title = "PD-associated activated fibroblast program", subtitle = stat_subtitle,
       x = NULL, y = "Activated module score", color = NULL, shape = "Dataset") +
  panel_theme + theme(legend.position = "bottom")

# E/F: mechanism-localization panels. Show expression-positive cells explicitly and retain zeros.
umap <- as.data.frame(Embeddings(fibro, "umap_fibro_final"))
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
expr <- FetchData(fibro, vars = c("AGT", "C3"), layer = "data")
plot_data <- cbind(umap, fibro[[]][rownames(umap), c("state_final", "group", "sample")], expr[rownames(umap), ])

pi16_n <- sum(plot_data$state_final == "PI16_Fib")
agt_n <- sum(plot_data$state_final == "PI16_Fib" & plot_data$AGT > 0)
agt_samples <- dplyr::n_distinct(plot_data$sample[plot_data$state_final == "PI16_Fib" & plot_data$AGT > 0])
p_e <- ggplot(plot_data, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey90", size = 0.08) +
  geom_point(data = filter(plot_data, state_final == "PI16_Fib"),
             color = state_colors[["PI16_Fib"]], size = 0.18, alpha = 0.65) +
  geom_point(data = filter(plot_data, state_final == "PI16_Fib", AGT > 0),
             aes(color = AGT), size = 1.4) +
  scale_color_viridis_c(option = "magma", begin = 0.35, end = 1) +
  labs(title = "AGT expression in PI16 fibroblasts",
       subtitle = sprintf("%d/%d cells across %d samples; descriptive localization", agt_n, pi16_n, agt_samples),
       x = "UMAP 1", y = "UMAP 2", color = "AGT") +
  coord_equal() + panel_theme

act_n <- sum(plot_data$state_final == "Activated_Fib")
c3_n <- sum(plot_data$state_final == "Activated_Fib" & plot_data$C3 > 0)
c3_samples <- dplyr::n_distinct(plot_data$sample[plot_data$state_final == "Activated_Fib" & plot_data$C3 > 0])
p_f <- ggplot(plot_data, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey90", size = 0.08) +
  geom_point(data = filter(plot_data, state_final == "Activated_Fib"),
             color = state_colors[["Activated_Fib"]], size = 0.18, alpha = 0.45) +
  geom_point(data = filter(plot_data, state_final == "Activated_Fib", C3 > 0),
             aes(color = C3), size = 0.28, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 1) +
  labs(title = "C3 expression in activated fibroblasts",
       subtitle = sprintf("%d/%d cells across %d samples", c3_n, act_n, c3_samples),
       x = "UMAP 1", y = "UMAP 2", color = "C3") +
  coord_equal() + panel_theme

save_plot(file.path(main_dir, "F6A_final_fibroblast_state_umap.pdf"), p_a, 8.5, 6.5)
save_plot(file.path(main_dir, "F6B_final_fixed_marker_dotplot.pdf"), p_b, 13, 5.6)
save_plot(file.path(main_dir, "F6C_final_sample_level_composition.pdf"), p_c, 17, 6.4)
save_plot(file.path(main_dir, "F6D_final_activated_program.pdf"), p_d, 7.2, 5.7)
save_plot(file.path(main_dir, "F6E_final_AGT_in_PI16.pdf"), p_e, 7.0, 5.7)
save_plot(file.path(main_dir, "F6F_final_C3_in_activated.pdf"), p_f, 7.0, 5.7)

# Journal-ready composite. Composition gets a full-width row because sample labels are dense.
figure6 <- ((p_a + p_b) + plot_layout(widths = c(0.85, 1.35))) /
  p_c /
  (p_d + p_e + p_f) +
  plot_layout(heights = c(1.0, 0.85, 0.9)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 20)))
save_plot(file.path(main_dir, "Figure6_fibroblast_revision_final.pdf"), figure6, 20, 18)

# Supplement: program-versus-stress specificity at the sample level.
supp_program <- sample_program |>
  select(sample, dataset, group, score_Activated_fibroblast, score_Stress, score_Contamination) |>
  pivot_longer(starts_with("score_"), names_to = "program", values_to = "score") |>
  mutate(program = recode(program,
                          score_Activated_fibroblast = "Activated program",
                          score_Stress = "Stress program",
                          score_Contamination = "Contamination program"))
p_specificity <- ggplot(supp_program, aes(group, score, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.12) +
  geom_point(aes(shape = dataset), position = position_jitter(width = 0.10), size = 1.9) +
  facet_wrap(~ program, scales = "free_y", nrow = 1) +
  scale_color_manual(values = group_colors) +
  labs(title = "Activated, stress, and contamination programs at sample level",
       subtitle = "The contamination program was not different after dataset adjustment (FDR=0.529)",
       x = NULL, y = "Mean module score", color = NULL, shape = "Dataset") +
  panel_theme + theme(legend.position = "bottom")
save_plot(file.path(supp_dir, "S08_program_stress_contamination_specificity.pdf"), p_specificity, 12, 4.8)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_figure_assembly.txt"))
message("Final Figure 6 assembled: ", file.path(main_dir, "Figure6_fibroblast_revision_final.pdf"))
