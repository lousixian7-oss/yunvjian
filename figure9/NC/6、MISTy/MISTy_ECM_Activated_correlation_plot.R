#!/usr/bin/env Rscript

.libPaths(c("G:/gurobi/Rlib45", .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

weights_file <- "RCTD_refinedFib_ECM_weights_full.csv"
output_pdf <- "ECM_Fib_Activated_Fib_spot_correlation.pdf"
output_csv <- "ECM_Fib_Activated_Fib_spot_correlation_stats.csv"

stopifnot(file.exists(weights_file))
d <- read.csv(weights_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("spot_id", "ECM_Fib", "Activated_Fib") %in% colnames(d)))

keep <- is.finite(d$ECM_Fib) & is.finite(d$Activated_Fib)
plot_df <- d[keep, c("spot_id", "ECM_Fib", "Activated_Fib")]

pearson <- cor.test(plot_df$ECM_Fib, plot_df$Activated_Fib, method = "pearson")
spearman <- suppressWarnings(
  cor.test(plot_df$ECM_Fib, plot_df$Activated_Fib, method = "spearman", exact = FALSE)
)

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 2.2e-16) return("< 2.2e-16")
  format.pval(p, digits = 3, eps = 1e-16)
}

stats_out <- data.frame(
  n_spots = nrow(plot_df),
  pearson_r = unname(pearson$estimate),
  pearson_p = pearson$p.value,
  spearman_rho = unname(spearman$estimate),
  spearman_p = spearman$p.value
)
write.csv(stats_out, output_csv, row.names = FALSE)

stat_label <- sprintf(
  "n = %d\nPearson r = %.3f, p %s\nSpearman rho = %.3f, p %s",
  nrow(plot_df),
  unname(pearson$estimate),
  format_p(pearson$p.value),
  unname(spearman$estimate),
  format_p(spearman$p.value)
)

p <- ggplot(plot_df, aes(x = ECM_Fib, y = Activated_Fib)) +
  geom_point(size = 1.35, alpha = 0.48, color = "#2878B5") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9,
              color = "#C82423", fill = "#F3B6B2", alpha = 0.28) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = stat_label,
    hjust = 1.06,
    vjust = 1.12,
    size = 3.8,
    linewidth = 0.25,
    fill = "white"
  ) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  labs(
    title = "ECM_Fib vs Activated_Fib across spatial spots",
    subtitle = "RCTD weights_full; each point represents one Visium spot",
    x = "ECM_Fib proportion",
    y = "Activated_Fib proportion",
    caption = "Correlation is spot-level and is not adjusted for spatial autocorrelation."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#555555"),
    plot.caption = element_text(color = "#666666", hjust = 0),
    axis.title = element_text(face = "bold")
  )

ggsave(output_pdf, p, width = 7.2, height = 5.8, device = cairo_pdf)
cat("Saved:", output_pdf, "\n")
print(stats_out)
