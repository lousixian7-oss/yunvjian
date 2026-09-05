args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1]] else
  "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final"
output_dir <- if (length(args) >= 2L) args[[2]] else
  "G:/1Yunvjian/FigureS4_fibroblast_dataset_bias"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(project_dir, "results", "tables")

state_levels <- c("Activated_Fib", "ECM_Fib", "PI16_Fib",
                  "Inflammatory_Fib", "Adventitial_Fib")
dataset_levels <- c("GSE152042", "GSE164241", "GSE171213")

# -----------------------------------------------------------------------------
# 1. Stress-adjusted activated-program sensitivity models (sample is the unit).
# -----------------------------------------------------------------------------
program <- read.csv(file.path(table_dir, "17_sample_level_module_scores.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
eligible <- subset(program,
  total_fibro >= 50 & dataset %in% c("GSE164241", "GSE171213"))
eligible$group <- factor(eligible$group, levels = c("NC", "PD"))
eligible$dataset <- factor(eligible$dataset)

gse164 <- subset(eligible, dataset == "GSE164241")
gse164$cohort <- factor(gse164$cohort)

extract_pd <- function(model, label, scope) {
  sm <- summary(model)$coefficients
  ci <- confint(model, "groupPD", level = 0.95)
  data.frame(
    scope = scope,
    model = label,
    effect = unname(sm["groupPD", "Estimate"]),
    std_error = unname(sm["groupPD", "Std. Error"]),
    t_value = unname(sm["groupPD", "t value"]),
    p_value = unname(sm["groupPD", "Pr(>|t|)"]),
    lower_95 = unname(ci[1]),
    upper_95 = unname(ci[2]),
    n_samples = stats::nobs(model),
    formula = paste(deparse(formula(model)), collapse = ""),
    stringsAsFactors = FALSE
  )
}

model_table <- rbind(
  extract_pd(
    lm(score_Activated_fibroblast ~ dataset + group, data = eligible),
    "Dataset adjusted", "GSE164241 + GSE171213"),
  extract_pd(
    lm(score_Activated_fibroblast ~ dataset + score_Stress + group,
       data = eligible),
    "Dataset + stress adjusted", "GSE164241 + GSE171213"),
  extract_pd(
    lm(score_Activated_fibroblast ~ cohort + group, data = gse164),
    "Cohort adjusted", "GSE164241"),
  extract_pd(
    lm(score_Activated_fibroblast ~ cohort + score_Stress + group,
       data = gse164),
    "Cohort + stress adjusted", "GSE164241")
)
model_table$significance <- cut(
  model_table$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns"), right = FALSE
)
write.csv(model_table,
          file.path(output_dir, "S4G2_stress_adjusted_activated_models.csv"),
          row.names = FALSE)
write.csv(model_table,
          file.path(table_dir, "27_stress_adjusted_activated_models.csv"),
          row.names = FALSE)

draw_stress_plot <- function(device = c("pdf", "png")) {
  device <- match.arg(device)
  out <- file.path(output_dir, paste0(
    "S4G2_stress_adjusted_activated_model.", device))
  if (device == "pdf") {
    grDevices::cairo_pdf(out, width = 9.0, height = 5.5)
  } else {
    grDevices::png(out, width = 2700, height = 1650, res = 300)
  }
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  par(mar = c(5.2, 10.5, 4.2, 1.5), family = "sans", las = 1)
  y <- rev(seq_len(nrow(model_table)))
  xlim <- range(c(model_table$lower_95, model_table$upper_95, 0), na.rm = TRUE)
  xpad <- diff(xlim) * 0.10
  plot_xlim <- c(xlim[1] - xpad, xlim[2] + max(0.24, xpad * 2.5))
  plot(model_table$effect, y, type = "n",
       xlim = plot_xlim, ylim = c(0.5, nrow(model_table) + 0.5),
       yaxt = "n", ylab = "", xlab = "PD effect on activated module score (95% CI)",
       main = "Stress-adjusted activated fibroblast program sensitivity",
       cex.main = 1.30, cex.lab = 1.15, cex.axis = 1.0)
  abline(v = 0, lty = 2, col = "grey55", lwd = 1.2)
  axis(2, at = y,
       labels = paste0(model_table$scope, "\n", model_table$model),
       tick = FALSE, cex.axis = 0.88)
  cols <- ifelse(grepl("stress", model_table$model, ignore.case = TRUE),
                 "#B2182B", "#3C78D8")
  segments(model_table$lower_95, y, model_table$upper_95, y,
           lwd = 2.3, col = cols)
  points(model_table$effect, y, pch = 19, cex = 1.35, col = cols)
  text(model_table$upper_95 + xpad * 0.15, y,
       labels = paste0(model_table$significance, "  P=",
         format.pval(model_table$p_value, digits = 2, eps = 1e-4)),
       pos = 4, cex = 0.92, font = 2)
  mtext("Overall effect attenuated after stress adjustment; the GSE164241 internal-cohort effect remained significant.",
        side = 3, line = 0.25, cex = 0.88, col = "grey30")
  legend("topright", legend = c("Without stress covariate", "Stress adjusted"),
         col = c("#3C78D8", "#B2182B"), pch = 19, lwd = 2,
         bty = "n", cex = 0.90)
}
draw_stress_plot("pdf")
draw_stress_plot("png")

# Revised Figure 6I: replace pooled/mean state composition with a continuous
# sample-level activated-program display. Low-cell samples remain display-only.
figure6_dir <- file.path(project_dir, "results", "figures",
                         "figure6_revised_GK")
dir.create(figure6_dir, recursive = TRUE, showWarnings = FALSE)
program$group <- factor(program$group, levels = c("NC", "PD"))
within_dataset_stats <- do.call(rbind, lapply(dataset_levels, function(ds) {
  z <- subset(program, dataset == ds & total_fibro >= 50)
  p <- if (length(unique(z$group)) == 2L &&
           all(table(z$group) >= 2L)) {
    wilcox.test(score_Activated_fibroblast ~ group, data = z,
                exact = FALSE)$p.value
  } else NA_real_
  data.frame(dataset = ds, n_eligible = nrow(z), p_value = p,
             stringsAsFactors = FALSE)
}))
write.csv(within_dataset_stats,
          file.path(figure6_dir,
                    "Figure6I_within_dataset_Wilcoxon_statistics.csv"),
          row.names = FALSE)

draw_figure6i <- function(device = c("pdf", "png")) {
  device <- match.arg(device)
  out <- file.path(figure6_dir,
    paste0("Figure6I_sample_level_activated_program.", device))
  if (device == "pdf") {
    grDevices::cairo_pdf(out, width = 7.2, height = 7.2)
  } else {
    grDevices::png(out, width = 2160, height = 2160, res = 300)
  }
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  par(mar = c(5.0, 5.4, 5.0, 1.4), family = "sans")
  set.seed(20260823)
  yr <- range(program$score_Activated_fibroblast, na.rm = TRUE)
  ypad <- max(diff(yr) * 0.25, 0.12)
  plot(NA, xlim = c(0.55, 2.45), ylim = c(yr[1] - 0.04, yr[2] + ypad),
       xaxt = "n", xlab = "Disease group",
       ylab = "Sample-level activated module score",
       main = "PD-associated activated fibroblast program",
       cex.main = 1.28, cex.lab = 1.08, cex.axis = 0.96)
  axis(1, at = 1:2, labels = c("NC", "PD"), cex.axis = 1.0)
  vals <- split(eligible$score_Activated_fibroblast, eligible$group)
  vals <- vals[c("NC", "PD")]
  boxplot(vals, at = 1:2, add = TRUE, axes = FALSE, outline = FALSE,
          boxwex = 0.56, border = c("#43B5C1", "#E56B5D"),
          col = c("#43B5C133", "#E56B5D33"), lwd = 1.8)
  ds_cols <- c(GSE152042 = "#D55E00", GSE164241 = "#009E73",
               GSE171213 = "#3C78D8")
  for (i in seq_len(nrow(program))) {
    x0 <- if (program$group[i] == "NC") 1 else 2
    x0 <- x0 + runif(1, -0.12, 0.12)
    is_eligible <- program$total_fibro[i] >= 50
    points(x0, program$score_Activated_fibroblast[i],
           pch = if (is_eligible) 19 else 1,
           cex = if (is_eligible) 1.08 else 1.18,
           lwd = 1.4, col = ds_cols[program$dataset[i]])
  }
  bracket_y <- yr[2] + ypad * 0.45
  segments(1, bracket_y, 2, bracket_y, lwd = 1.4)
  segments(1, bracket_y, 1, bracket_y - ypad * 0.08, lwd = 1.4)
  segments(2, bracket_y, 2, bracket_y - ypad * 0.08, lwd = 1.4)
  text(1.5, bracket_y + ypad * 0.08,
       "Primary: **; stress-adjusted: ns", cex = 0.82, font = 2)
  mtext("Each GSM is one replicate", side = 3, line = 1.15,
        cex = 0.80, col = "grey30")
  mtext("Primary dataset-adjusted FDR=0.00331; stress-adjusted P=0.103 (Fig. S4G)",
        side = 3, line = 0.20, cex = 0.72, col = "grey30")
  legend("topleft", legend = names(ds_cols), col = ds_cols, pch = 19,
         title = "Dataset", bty = "n", cex = 0.86, pt.cex = 1.0)
  legend("topright", legend = c(">=50 fibroblasts (model)",
                                "<50 fibroblasts (display only)"),
         pch = c(19, 1), col = "grey25", bty = "n", cex = 0.82)
  box(col = "grey35")
}
draw_figure6i("pdf")
draw_figure6i("png")

# -----------------------------------------------------------------------------
# 2. Resolution agreement and donor/GSM representation.
# -----------------------------------------------------------------------------
ari <- read.csv(file.path(table_dir, "05_resolution_stability_ARI.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
ari$res_a <- sub(".*res\\.", "", ari$a)
ari$res_b <- sub(".*res\\.", "", ari$b)
res_levels <- c("0.2", "0.3", "0.4", "0.5")
ari_mat <- matrix(NA_real_, nrow = 4, ncol = 4,
                  dimnames = list(res_levels, res_levels))
for (i in seq_len(nrow(ari))) {
  ari_mat[ari$res_b[i], ari$res_a[i]] <- ari$ARI[i]
}

composition <- read.csv(file.path(table_dir,
  "16_sample_level_fibroblast_composition.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
composition <- subset(composition, state %in% state_levels)
composition$eligible <- composition$total_fibro >= 50

support_rows <- list()
k <- 1L
for (ds in dataset_levels) {
  ds_dat <- subset(composition, dataset == ds & eligible)
  n_total <- length(unique(ds_dat$sample))
  for (st in state_levels) {
    st_dat <- subset(ds_dat, state == st)
    n_supported <- if (n_total == 0L) 0L else sum(st_dat$n >= 10)
    support_rows[[k]] <- data.frame(
      dataset = ds, state = st, supported_GSM = n_supported,
      eligible_GSM = n_total,
      support_fraction = if (n_total == 0L) NA_real_ else n_supported / n_total,
      label = if (n_total == 0L) "display only" else
        paste0(n_supported, "/", n_total),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
support <- do.call(rbind, support_rows)
write.csv(support,
          file.path(output_dir, "S4H_state_GSM_representation.csv"),
          row.names = FALSE)
write.csv(support,
          file.path(table_dir, "28_state_GSM_representation.csv"),
          row.names = FALSE)

palette_fun <- colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))
heat_cols <- palette_fun(101)

draw_tile_matrix <- function(values, row_labels, col_labels, text_labels,
                             zlim, title, subtitle, na_col = "grey90") {
  nr <- nrow(values); nc <- ncol(values)
  plot.new()
  plot.window(xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
              xaxs = "i", yaxs = "i")
  for (r in seq_len(nr)) for (c in seq_len(nc)) {
    val <- values[r, c]
    col <- if (is.na(val)) na_col else {
      idx <- round((val - zlim[1]) / diff(zlim) * 100) + 1
      heat_cols[max(1, min(101, idx))]
    }
    rect(c - 0.5, nr - r + 0.5, c + 0.5, nr - r + 1.5,
         col = col, border = "white", lwd = 1.5)
    text(c, nr - r + 1, labels = text_labels[r, c],
         cex = 0.92, font = 2,
         col = if (!is.na(val) && val > mean(zlim)) "white" else "black")
  }
  axis(1, at = seq_len(nc), labels = col_labels, tick = FALSE,
       cex.axis = 0.90)
  axis(2, at = seq_len(nr), labels = rev(row_labels), tick = FALSE,
       las = 1, cex.axis = 0.88)
  box(col = "grey40")
  title(main = title, cex.main = 1.18, font.main = 2)
  mtext(subtitle, side = 3, line = 0.25, cex = 0.78, col = "grey30")
}

draw_robustness <- function(device = c("pdf", "png")) {
  device <- match.arg(device)
  out <- file.path(output_dir, paste0("S4H_cluster_robustness.", device))
  if (device == "pdf") {
    grDevices::cairo_pdf(out, width = 12.0, height = 5.8)
  } else {
    grDevices::png(out, width = 3600, height = 1740, res = 300)
  }
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  layout(matrix(1:2, nrow = 1), widths = c(0.88, 1.25))

  par(mar = c(4.6, 4.8, 4.2, 1.0), family = "sans")
  draw_tile_matrix(
    ari_mat, res_levels, res_levels,
    matrix(sprintf("%.2f", ari_mat), 4, 4),
    zlim = c(0.5, 1.0),
    title = "Partition stability",
    subtitle = "Adjusted Rand index across resolution 0.2-0.5"
  )

  support_mat <- matrix(NA_real_, nrow = length(state_levels),
                        ncol = length(dataset_levels),
                        dimnames = list(state_levels, dataset_levels))
  label_mat <- matrix("", nrow = length(state_levels),
                      ncol = length(dataset_levels),
                      dimnames = list(state_levels, dataset_levels))
  for (i in seq_len(nrow(support))) {
    support_mat[support$state[i], support$dataset[i]] <-
      support$support_fraction[i]
    label_mat[support$state[i], support$dataset[i]] <- support$label[i]
  }
  par(mar = c(4.6, 8.0, 4.2, 1.0), family = "sans")
  draw_tile_matrix(
    support_mat, state_levels, dataset_levels, label_mat,
    zlim = c(0, 1),
    title = "Donor representation",
    subtitle = "GSMs with >=10 state cells / eligible GSMs; GSE152042 is display-only"
  )
}
draw_robustness("pdf")
draw_robustness("png")

capture.output(sessionInfo(),
  file = file.path(output_dir, "sessionInfo_mandatory_validations.txt"))
message("Mandatory validation outputs written to: ", output_dir)
