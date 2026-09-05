args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1L) args[[1]] else "G:/1Yunvjian/0a26.7.7singlecell/8.10"
output_dir <- if (length(args) >= 2L) args[[2]] else file.path(root_dir, "fibroblast_revision_final")

source_dir <- file.path(output_dir, "results/figures")
supp_dir <- file.path(source_dir, "reviewer_dataset_bias_supplement")
log_dir <- file.path(output_dir, "logs")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(log_dir, "06_assemble_reviewer_supplement.log"), "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_con)
}, add = TRUE)

suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
})

image_panel <- function(path, tag) {
  stopifnot(file.exists(path))
  img <- png::readPNG(path)
  grob <- grobTree(
    rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE),
    textGrob(tag, x = unit(0.008, "npc"), y = unit(0.992, "npc"),
             just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 24))
  )
  wrap_elements(full = grob)
}

p_a <- image_panel(file.path(source_dir, "supplement/S01_contamination_marker_dotplot.png"), "A")
p_b <- image_panel(file.path(source_dir, "supplement/S02_fibroblast_cleaning_before_after_umap.png"), "B")

fib_colors <- c(
  ECM_Fib = "#F8766D", PI16_Fib = "#CCCC00",
  Inflammatory_Fib = "#00BFC4", Activated_Fib = "#3399FF",
  Adventitial_Fib = "#C77CFF"
)
composition <- read.csv(
  file.path(output_dir,
            "results/tables/figure6_revised_GK/I_all_GSM_sample_state_composition.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
) %>%
  mutate(
    state = factor(state, levels = names(fib_colors)),
    sample_label = paste0(sample, "\n(n=", total_fibro, ")")
  )

# C: all biological samples are displayed, including low-cell GSE152042 samples.
p_c <- ggplot(composition,
              aes(x = sample_label, y = fraction, fill = state)) +
  geom_col(width = 0.86, color = "white", linewidth = 0.15) +
  facet_grid(. ~ dataset + group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fib_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.01))) +
  labs(
    tag = "C", x = NULL, y = "Fibroblast-state fraction", fill = "state",
    title = "Composition by biological sample (all three GSEs)",
    subtitle = "One bar per GSM; n is the qualified-fibroblast denominator"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    axis.text.x = element_text(angle = 70, hjust = 1, size = 7),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right"
  )

# D: reassigned five-state object; no unresolved state is shown.
contingency <- xtabs(n_state ~ dataset + state, data = composition)
chi <- suppressWarnings(chisq.test(contingency, correct = FALSE))
cramers_v <- sqrt(as.numeric(chi$statistic) /
                    (sum(contingency) * min(nrow(contingency) - 1,
                                            ncol(contingency) - 1)))
heat_df <- as.data.frame(prop.table(contingency, margin = 2))
colnames(heat_df) <- c("dataset", "state", "within_state_fraction")
p_d <- ggplot(heat_df, aes(x = dataset, y = state, fill = within_state_fraction)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = scales::percent(within_state_fraction, accuracy = 0.1)),
            size = 3.8) +
  scale_fill_viridis_c(labels = scales::percent, limits = c(0, 1)) +
  labs(
    tag = "D", x = NULL, y = NULL, fill = "Within-state\nfraction",
    title = "Dataset contribution to the final five fibroblast states",
    subtitle = sprintf("Cramer's V = %.3f; residual dataset association is retained in downstream models",
                       cramers_v)
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(supp_dir, "panelC_all_GSM_composition.png"), p_c,
       width = 16, height = 6.2, dpi = 300, limitsize = FALSE)
ggsave(file.path(supp_dir, "panelD_final_five_state_dataset_heatmap.png"), p_d,
       width = 8.5, height = 5.2, dpi = 300, limitsize = FALSE)
p_e <- image_panel(file.path(source_dir, "supplement/S04_activated_program_by_dataset.png"), "E")
p_f <- image_panel(file.path(source_dir, "supplement/S05_GSE164241_activated_program_cohort_sensitivity.png"), "F")
p_g <- image_panel(file.path(source_dir, "supplement/S08_program_stress_contamination_specificity.png"), "G")

reviewer_supp <- p_a / p_b / p_c / (p_d + p_e) / (p_f + p_g) +
  plot_layout(heights = c(0.72, 0.86, 1.02, 0.78, 0.78))

pdf_path <- file.path(supp_dir, "FigureS_fibroblast_dataset_bias_reviewer_response.pdf")
png_path <- file.path(supp_dir, "FigureS_fibroblast_dataset_bias_reviewer_response.png")
ggsave(pdf_path, reviewer_supp, width = 19, height = 22,
       device = cairo_pdf, limitsize = FALSE)
ggsave(png_path, reviewer_supp, width = 19, height = 22,
       dpi = 300, limitsize = FALSE)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_06_reviewer_supplement.txt"))
message("Completed reviewer-response supplement: ", pdf_path)
