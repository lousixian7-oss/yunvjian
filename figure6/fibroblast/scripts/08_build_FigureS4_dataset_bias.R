args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1]] else
  "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final"
output_dir <- if (length(args) >= 2L) args[[2]] else
  "G:/1Yunvjian/FigureS4_fibroblast_dataset_bias"

local_library <- file.path(dirname(project_dir),
  "batch_integration_revision_20260821", "R_library")
if (dir.exists(local_library)) {
  .libPaths(c(local_library, .libPaths()))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(limma)
  library(ggpubr)
})

object_file <- file.path(project_dir,
  "results/objects/fibroblast_revision_final_no_unassigned.rds")
table_dir <- file.path(project_dir, "results/tables")
figure6_table_dir <- file.path(table_dir, "figure6_revised_GK")

fibro <- readRDS(object_file)
DefaultAssay(fibro) <- "RNA"
fibro[["RNA"]] <- JoinLayers(fibro[["RNA"]])

state_levels <- c("Activated_Fib", "ECM_Fib", "PI16_Fib",
                  "Inflammatory_Fib", "Adventitial_Fib")
fibro$state_final <- factor(as.character(fibro$state_final), levels = state_levels)
dataset_colors <- c(GSE152042 = "#D55E00", GSE164241 = "#009E73",
                    GSE171213 = "#3C78D8")
group_colors <- c(NC = "#43B5C1", PD = "#E56B5D")
state_colors <- c(Activated_Fib = "#247FD1", ECM_Fib = "#E85C50",
                  PI16_Fib = "#B3B300", Inflammatory_Fib = "#00A5AE",
                  Adventitial_Fib = "#AE63DB")

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1,
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25",
                                   margin = margin(b = 6)),
      axis.title = element_text(face = "plain", size = base_size),
      axis.text = element_text(size = base_size - 1, color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey65",
                                      linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.position = "right",
      plot.margin = margin(7, 7, 7, 7)
    )
}

save_panel <- function(short_name, descriptive_name, plot, width, height) {
  short_pdf <- file.path(output_dir, paste0(short_name, ".pdf"))
  descriptive_pdf <- file.path(output_dir, paste0(descriptive_name, ".pdf"))
  ggsave(short_pdf, plot, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  file.copy(short_pdf, descriptive_pdf, overwrite = TRUE)
  ggsave(file.path(output_dir, paste0(short_name, ".png")), plot,
         width = width, height = height, dpi = 300, limitsize = FALSE)
}

add_panel_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag = element_text(face = "bold", size = 22),
          plot.tag.position = c(0.01, 0.99))
}

# S4A: fixed contamination-marker panel on the cleaned fibroblast object.
contamination_markers <- list(
  `Pericyte / smooth muscle` = c("RGS5", "GUCY1A1", "GUCY1B1", "MCAM"),
  `T / NK` = c("CD3D", "TRAC", "NKG7", "GZMK"),
  Epithelial = c("KRT5", "KRT6A", "KRT13", "DSP"),
  Plasma = c("MZB1", "CD79A", "IGHG", "IGKC"),
  Myeloid = c("CD68", "ACP5")
)
contamination_markers <- lapply(contamination_markers, intersect,
                                y = rownames(fibro))
contamination_markers <- contamination_markers[lengths(contamination_markers) > 0]
Idents(fibro) <- fibro$state_final
p_a <- DotPlot(fibro, features = contamination_markers, assay = "RNA",
               dot.scale = 6.5, cols = c("grey92", "#B2182B")) +
  RotatedAxis() +
  labs(x = NULL, y = NULL, color = "Average\nexpression",
       size = "Percent\nexpressed",
       title = "Fibroblast purity validation",
       subtitle = "Lineage-contamination markers after fibroblast cleaning") +
  theme_pub(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        panel.spacing.x = unit(0.25, "lines"))

# S4B: before/after integration, colored by dataset. Disease group is not used.
before_coords <- as.data.frame(Embeddings(fibro, "umap_fibro_before_final"))
before_coords$dataset <- as.character(fibro$dataset)
colnames(before_coords)[1:2] <- c("x", "y")
after_coords <- as.data.frame(Embeddings(fibro, "umap_fibro_final"))
after_coords$dataset <- as.character(fibro$dataset)
colnames(after_coords)[1:2] <- c("x", "y")
rare_before <- before_coords %>% filter(dataset == "GSE152042")
rare_after <- after_coords %>% filter(dataset == "GSE152042")
p_b1 <- DimPlot(fibro, reduction = "umap_fibro_before_final",
                group.by = "dataset", cols = dataset_colors,
                pt.size = 0.90, shuffle = TRUE, raster = TRUE,
                raster.dpi = c(600, 600)) +
  geom_point(data = rare_before, aes(x = x, y = y), inherit.aes = FALSE,
             shape = 21, fill = dataset_colors[["GSE152042"]], color = "black",
             size = 1.15, stroke = 0.22) +
  labs(title = "Before Harmony", x = "UMAP 1", y = "UMAP 2", color = "Dataset") +
  theme_pub(12) +
  theme(aspect.ratio = 1, panel.border = element_rect(fill = NA,
        color = "black", linewidth = 0.45), legend.position = "none")
p_b2 <- DimPlot(fibro, reduction = "umap_fibro_final",
                group.by = "dataset", cols = dataset_colors,
                pt.size = 0.90, shuffle = TRUE, raster = TRUE,
                raster.dpi = c(600, 600)) +
  geom_point(data = rare_after, aes(x = x, y = y), inherit.aes = FALSE,
             shape = 21, fill = dataset_colors[["GSE152042"]], color = "black",
             size = 1.15, stroke = 0.22) +
  labs(title = "After Harmony", subtitle = "Dataset corrected by Harmony",
       x = "UMAP 1", y = "UMAP 2", color = "Dataset") +
  theme_pub(12) +
  theme(aspect.ratio = 1, panel.border = element_rect(fill = NA,
        color = "black", linewidth = 0.45))
p_b <- p_b1 + p_b2 + plot_layout(guides = "collect") &
  theme(legend.position = "right")

# S4C: every GSM is one point; low-cell samples are display-only open symbols.
composition <- read.csv(file.path(figure6_table_dir,
  "I_all_GSM_sample_state_composition.csv"), stringsAsFactors = FALSE,
  check.names = FALSE) %>%
  mutate(
    state = factor(state, levels = state_levels),
    group = factor(group, levels = c("NC", "PD")),
    eligible = total_fibro >= 50
  )
# Refit the final five-state composition model after the rare audit state was
# assigned to ECM_Fib. Each GSM is one column; dataset is retained in the design.
eligible_composition <- composition %>% filter(eligible)
composition_sample_meta <- eligible_composition %>%
  distinct(sample, group, dataset) %>% arrange(sample) %>%
  mutate(group = droplevels(group), dataset = factor(dataset))
composition_counts <- xtabs(n_state ~ state + sample, data = eligible_composition)
composition_counts <- composition_counts[, composition_sample_meta$sample, drop = FALSE]
composition_prop <- sweep(composition_counts, 2, colSums(composition_counts), "/")
composition_design <- model.matrix(~ dataset + group, data = composition_sample_meta)
composition_fit <- eBayes(
  lmFit(asin(sqrt(composition_prop)), composition_design),
  robust = TRUE, trend = TRUE
)
composition_stats <- topTable(
  composition_fit, coef = "groupPD", number = Inf, sort.by = "none"
) %>%
  tibble::rownames_to_column("state") %>%
  transmute(state, effect = logFC, p_value = P.Value, FDR = adj.P.Val) %>%
  mutate(
    state = factor(state, levels = state_levels),
    group1 = "NC", group2 = "PD",
    significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  left_join(
    composition %>% group_by(state) %>%
      summarise(y.position = min(1.08, max(fraction, na.rm = TRUE) * 1.08 + 0.02),
                .groups = "drop"),
    by = "state"
  )
write.csv(composition_stats,
          file.path(output_dir, "S4C_dataset_adjusted_composition_statistics.csv"),
          row.names = FALSE)

p_c <- ggplot(composition, aes(x = group, y = fraction)) +
  geom_boxplot(data = filter(composition, eligible), aes(fill = group),
               width = 0.56, outlier.shape = NA, alpha = 0.20,
               color = "grey25", linewidth = 0.45) +
  geom_point(aes(color = dataset, shape = eligible),
             position = position_jitter(width = 0.11, height = 0),
             size = 2.55, stroke = 0.85) +
  ggpubr::stat_pvalue_manual(
    composition_stats, label = "significance",
    xmin = "group1", xmax = "group2", y.position = "y.position",
    inherit.aes = FALSE, tip.length = 0.012, bracket.size = 0.45,
    size = 4.2, hide.ns = FALSE
  ) +
  facet_wrap(~ state, nrow = 1, scales = "free_y") +
  scale_color_manual(values = dataset_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = ">=50 fibroblasts (model)",
               `FALSE` = "<50 fibroblasts (display only)")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = "Disease group", y = "Fraction of fibroblast state",
       color = "Dataset", fill = NULL, shape = "Sample inclusion",
       title = "Sample-level fibroblast-state proportions",
       subtitle = "All three datasets shown; stars denote dataset-adjusted FDR") +
  theme_pub(12) +
  theme(legend.position = "right")

# S4D: retain and report moderate residual dataset association.
contingency <- xtabs(n_state ~ dataset + state, data = composition)
chi <- suppressWarnings(chisq.test(contingency, correct = FALSE))
cramers_v <- sqrt(as.numeric(chi$statistic) /
                    (sum(contingency) * min(nrow(contingency) - 1,
                                            ncol(contingency) - 1)))
heat_df <- as.data.frame(prop.table(contingency, margin = 2))
colnames(heat_df) <- c("dataset", "state", "within_state_fraction")
p_d <- ggplot(heat_df,
              aes(x = dataset, y = state, fill = within_state_fraction)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = scales::percent(within_state_fraction, accuracy = 0.1)),
            size = 3.8, fontface = "bold") +
  scale_fill_viridis_c(option = "C", limits = c(0, 1),
                       labels = scales::percent) +
  labs(x = NULL, y = NULL, fill = "Within-state\ndataset fraction",
       title = "Dataset contribution and residual association",
       subtitle = paste0("Cramer's V = ", sprintf("%.3f", cramers_v),
         " (moderate); disease effects were tested after dataset adjustment")) +
  theme_pub(12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# S4E-G use sample-level means from uncorrected normalized expression scores.
sample_program <- read.csv(file.path(table_dir,
  "17_sample_level_module_scores.csv"), stringsAsFactors = FALSE,
  check.names = FALSE) %>%
  mutate(group = factor(group, levels = c("NC", "PD")),
         dataset = factor(dataset,
           levels = c("GSE152042", "GSE164241", "GSE171213")),
         eligible = total_fibro >= 50)
eligible_program <- sample_program %>%
  filter(eligible, dataset %in% c("GSE164241", "GSE171213"))

# Pre-specified sensitivity analysis: test whether the disease-associated
# activated program remains after accounting for the sample-level stress score.
# This is a sensitivity model, not a claim that stress is biologically irrelevant.
extract_group_effect <- function(formula, model_name) {
  design <- model.matrix(formula, data = eligible_program)
  y <- matrix(eligible_program$score_Activated_fibroblast, nrow = 1,
              dimnames = list("Activated", eligible_program$sample))
  fit <- eBayes(lmFit(y, design), robust = TRUE, trend = TRUE)
  coef_name <- "groupPD"
  coef_index <- match(coef_name, colnames(design))
  estimate <- fit$coefficients[1, coef_index]
  se <- fit$stdev.unscaled[1, coef_index] * sqrt(fit$s2.post[1])
  crit <- qt(0.975, df = fit$df.total[1])
  data.frame(
    model = model_name,
    effect = estimate,
    lower = estimate - crit * se,
    upper = estimate + crit * se,
    p_value = fit$p.value[1, coef_index],
    n_samples = ncol(y),
    formula = paste(deparse(formula), collapse = "")
  )
}

stress_adjusted_stats <- bind_rows(
  extract_group_effect(~ dataset + group, "Dataset-adjusted"),
  extract_group_effect(~ dataset + score_Stress + group,
                       "Dataset + stress adjusted")
) %>%
  mutate(
    model = factor(model,
      levels = c("Dataset-adjusted", "Dataset + stress adjusted")),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
write.csv(stress_adjusted_stats,
          file.path(output_dir, "S4G_stress_adjusted_activated_model.csv"),
          row.names = FALSE)

p_e <- ggplot(sample_program,
              aes(x = group, y = score_Activated_fibroblast, color = group)) +
  geom_boxplot(data = eligible_program, outlier.shape = NA,
               width = 0.56, alpha = 0.15,
               linewidth = 0.5) +
  geom_point(aes(shape = eligible), position = position_jitter(width = 0.10),
             size = 2.7, alpha = 0.95, stroke = 0.85) +
  ggpubr::stat_compare_means(
    data = eligible_program, method = "wilcox.test", label = "p.signif",
    hide.ns = FALSE, size = 4.2, label.y.npc = 0.94
  ) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_y") +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = ">=50 fibroblasts (model)",
               `FALSE` = "<50 fibroblasts (display only)")) +
  labs(x = "Disease group", y = "Sample-level activated module score",
       color = NULL, shape = "Sample inclusion",
       title = "Dataset-adjusted activated fibroblast validation",
       subtitle = paste0(
         "All 3 datasets shown; GSE152042 display-only. ",
         "Overall adjusted FDR = 0.00331 (**)"
       )) +
  theme_pub(12) + theme(legend.position = "right")

gse164 <- eligible_program %>% filter(dataset == "GSE164241")
p_f <- ggplot(gse164,
              aes(x = group, y = score_Activated_fibroblast,
                  color = group, shape = cohort)) +
  geom_boxplot(aes(group = group), outlier.shape = NA, width = 0.56,
               alpha = 0.15, linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.10), size = 2.8,
             alpha = 0.95) +
  ggpubr::stat_pvalue_manual(
    data.frame(group1 = "NC", group2 = "PD",
               y.position = max(gse164$score_Activated_fibroblast) * 1.08,
               significance = "**"),
    label = "significance", xmin = "group1", xmax = "group2",
    y.position = "y.position", inherit.aes = FALSE,
    tip.length = 0.012, bracket.size = 0.45, size = 4.5
  ) +
  scale_color_manual(values = group_colors) +
  labs(x = "Disease group", y = "Sample-level activated module score",
       color = NULL, shape = "Internal cohort",
       title = "GSE164241 internal-cohort sensitivity analysis",
       subtitle = "Cohort-adjusted effect size = 0.357; FDR = 0.00257 (**)") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.16))) +
  theme_pub(12)

program_long <- eligible_program %>%
  select(sample, group, dataset,
         Activated = score_Activated_fibroblast,
         Stress = score_Stress,
         Contamination = score_Contamination) %>%
  pivot_longer(c(Activated, Stress, Contamination),
               names_to = "program", values_to = "score") %>%
  mutate(program = factor(program,
    levels = c("Activated", "Stress", "Contamination")))
program_labs <- c(
  Activated = "Activated\nFDR = 0.00331 (**)",
  Stress = "Stress\nFDR = 4.65e-6 (***)",
  Contamination = "Contamination\nFDR = 0.529 (ns)"
)
program_pvalues <- program_long %>%
  group_by(program) %>%
  summarise(y.position = max(score, na.rm = TRUE) * 1.08, .groups = "drop") %>%
  mutate(
    group1 = "NC", group2 = "PD",
    significance = c("**", "***", "ns")
  )
p_g_scores <- ggplot(program_long, aes(x = group, y = score, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.56, alpha = 0.15,
               linewidth = 0.5) +
  geom_point(aes(shape = dataset), position = position_jitter(width = 0.10),
             size = 2.45, alpha = 0.95) +
  ggpubr::stat_pvalue_manual(
    program_pvalues, label = "significance",
    xmin = "group1", xmax = "group2", y.position = "y.position",
    inherit.aes = FALSE, tip.length = 0.012, bracket.size = 0.45,
    size = 4.2, hide.ns = FALSE
  ) +
  facet_wrap(~ program, scales = "free_y", labeller = as_labeller(program_labs)) +
  scale_color_manual(values = group_colors) +
  labs(x = "Disease group", y = "Sample-level module score",
       color = NULL, shape = "Dataset",
       title = "Activated, stress and contamination program dissection",
       subtitle = "Activated and stress increased in PD; contamination did not differ") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.16))) +
  theme_pub(12)

stress_p <- stress_adjusted_stats %>%
  filter(model == "Dataset + stress adjusted") %>% pull(p_value)
p_g_adjusted <- ggplot(stress_adjusted_stats,
                       aes(x = effect, y = model)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55",
             linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.14,
                 linewidth = 0.65, color = "grey25") +
  geom_point(size = 3.2, color = "#B2182B") +
  geom_text(aes(label = significance), nudge_y = 0.18,
            size = 4.2, fontface = "bold") +
  labs(x = "PD effect on activated module score (95% CI)", y = NULL,
       title = "Activated program after stress adjustment",
       subtitle = paste0("Dataset + stress adjusted P = ",
                         format.pval(stress_p, digits = 3, eps = 1e-4))) +
  theme_pub(12) +
  theme(legend.position = "none")

p_g <- p_g_scores / p_g_adjusted +
  plot_layout(heights = c(1.7, 1.0))

# S4H: minimum cluster-robustness audit. ARI quantifies agreement between
# nearby graph partitions; donor support reports how many eligible GSM samples
# contribute at least 10 cells to each descriptive fibroblast state.
ari_long <- read.csv(file.path(table_dir, "05_resolution_stability_ARI.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    resolution_a = sub(".*res\\.", "", a),
    resolution_b = sub(".*res\\.", "", b),
    resolution_a = factor(resolution_a, levels = c("0.2", "0.3", "0.4", "0.5")),
    resolution_b = factor(resolution_b, levels = c("0.5", "0.4", "0.3", "0.2"))
  )
p_h1 <- ggplot(ari_long, aes(x = resolution_a, y = resolution_b, fill = ARI)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", ARI)), size = 3.6,
            fontface = "bold") +
  scale_fill_viridis_c(option = "C", limits = c(0.5, 1.0)) +
  labs(x = "Resolution", y = "Resolution", fill = "ARI",
       title = "Partition stability",
       subtitle = "Adjusted Rand index across resolution 0.2-0.5") +
  theme_pub(11)

donor_support <- composition %>%
  filter(eligible) %>%
  group_by(dataset, state) %>%
  summarise(
    supported_GSM = sum(n_state >= 10),
    eligible_GSM = n_distinct(sample),
    support_fraction = supported_GSM / eligible_GSM,
    .groups = "drop"
  ) %>%
  complete(
    dataset = factor(c("GSE152042", "GSE164241", "GSE171213"),
                     levels = c("GSE152042", "GSE164241", "GSE171213")),
    state = factor(state_levels, levels = state_levels),
    fill = list(supported_GSM = 0, eligible_GSM = 0,
                support_fraction = NA_real_)
  ) %>%
  mutate(label = ifelse(eligible_GSM == 0, "display only",
                        paste0(supported_GSM, "/", eligible_GSM)))
write.csv(donor_support,
          file.path(output_dir, "S4H_state_GSM_representation.csv"),
          row.names = FALSE)

p_h2 <- ggplot(donor_support,
               aes(x = dataset, y = state, fill = support_fraction)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = label), size = 3.25, fontface = "bold") +
  scale_fill_viridis_c(option = "D", limits = c(0, 1),
                       na.value = "grey90", labels = scales::percent) +
  labs(x = NULL, y = NULL, fill = "GSM support",
       title = "Donor representation",
       subtitle = "GSMs with >=10 state cells / eligible GSMs") +
  theme_pub(11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_h <- p_h1 | p_h2

save_panel("S4A", "S4A_fibroblast_purity", p_a, 10.2, 5.4)
save_panel("S4B", "S4B_fibroblast_Harmony", p_b, 11.2, 5.0)
save_panel("S4C", "S4C_sample_level_fraction", p_c, 14.2, 5.3)
save_panel("S4D", "S4D_dataset_association", p_d, 8.8, 5.5)
save_panel("S4E", "S4E_Activated_module_validation", p_e, 8.8, 5.2)
save_panel("S4F", "S4F_GSE164241_sensitivity", p_f, 7.4, 5.2)
save_panel("S4G", "S4G_program_dissection", p_g, 10.4, 8.0)
save_panel("S4H", "S4H_cluster_robustness", p_h, 11.2, 5.4)

p_a_comb <- p_a + theme(legend.position = "none",
                        strip.text = element_text(size = 7.5))
p_b_comb <- (add_panel_tag(p_b1 + theme(legend.position = "none"), "B") |
               (p_b2 +
                  theme(legend.position = "bottom", legend.direction = "horizontal",
                        legend.text = element_text(size = 6.5),
                        legend.title = element_text(size = 7.0)) +
                  guides(color = guide_legend(nrow = 1,
                    override.aes = list(size = 3)))))
p_c_comb <- p_c +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "vertical", legend.text = element_text(size = 6.2),
        legend.title = element_text(size = 6.8),
        strip.text = element_text(size = 7.5)) +
  guides(fill = "none", color = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1))
p_d_comb <- p_d + theme(legend.position = "none")
p_e_comb <- p_e + theme(legend.position = "none")
p_f_comb <- p_f + theme(legend.position = "none")
p_g_comb <- p_g + theme(legend.position = "none")
p_h_comb <- p_h + theme(legend.position = "none")

combined <-
  (add_panel_tag(p_a_comb, "A") | p_b_comb |
     add_panel_tag(p_c_comb, "C") | add_panel_tag(p_d_comb, "D")) /
  (add_panel_tag(p_e_comb, "E") | add_panel_tag(p_f_comb, "F") |
     add_panel_tag(p_g_comb, "G") | add_panel_tag(p_h_comb, "H")) +
  plot_layout(widths = c(1.25, 1.15, 1.35, 1.05),
              heights = c(1.0, 0.95))

ggsave(file.path(output_dir, "FigureS4_combined.pdf"), combined,
       width = 28.0, height = 13.6, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(output_dir, "FigureS4_combined.png"), combined,
       width = 28.0, height = 13.6, dpi = 300, limitsize = FALSE)

write.csv(data.frame(
  metric = c("Cramer's V", "Activated module FDR",
             "GSE164241 cohort-adjusted effect", "GSE164241 sensitivity FDR",
             "Stress FDR", "Contamination FDR",
             "Stress-adjusted activated P"),
  value = c(cramers_v, 0.00331011570134325, 0.357493209510718,
            0.00257376349286528, 4.65167662635933e-06,
            0.529017310885229, stress_p)
), file.path(output_dir, "FigureS4_key_statistics.csv"), row.names = FALSE)

writeLines(c(
  "Figure S4 analysis structure",
  "A. Tests whether the cleaned fibroblast object retains lineage-contamination signatures.",
  "B. Demonstrates de novo fibroblast re-clustering and dataset-only Harmony correction; disease group was not used.",
  "C. Replaces pooled-cell percentages with one observation per GSM. Samples with <50 fibroblasts are open symbols and display-only.",
  "D. Quantifies residual dataset-state association and retains Cramer's V rather than claiming complete batch removal.",
  "E. Tests the activated program across datasets using sample-level scores and a dataset-adjusted model.",
  "F. Tests whether the GSE164241 result survives adjustment for the GSM500/GSM517 internal cohort.",
  "G. Displays activated, stress and contamination scores and tests whether the activated disease effect persists after sample-level stress adjustment.",
  "H. Audits partition agreement across resolutions and the number of eligible GSM samples supporting each descriptive fibroblast state.",
  "Interpretation: PD-associated activated fibroblast programs remain after cleaning, sample-level quantification and dataset adjustment. The figure does not claim that Activated_Fib is the only PD-associated state or prove a causal transition."
), file.path(output_dir, "FigureS4_analysis_logic.txt"))

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_FigureS4.txt"))
message("Figure S4 completed in: ", output_dir)
