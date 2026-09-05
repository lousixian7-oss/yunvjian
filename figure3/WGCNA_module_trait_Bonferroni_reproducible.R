## Reproducible multiple-testing correction for the WGCNA module-trait panel
## Input:  .RData saved by the original WGCNA analysis
## Output: corrected results table and replacement Figure 3C panel

user_library <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
.libPaths(c(user_library, .libPaths()))

options(stringsAsFactors = FALSE)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

## On some Windows installations, commandArgs() may garble non-ASCII paths.
## Prefer the working directory whenever it contains the saved analysis state.
running_from_analysis_dir <- file.exists(".RData")
analysis_dir <- if (running_from_analysis_dir) {
  "."
} else {
  get_script_dir()
}
## Keep paths relative when running inside a non-ASCII Windows directory.
## R 4.5.x may fail when the same path is expanded to its full Chinese form.
input_file <- if (running_from_analysis_dir) ".RData" else file.path(analysis_dir, ".RData")
output_dir <- if (running_from_analysis_dir) "Bonferroni_correction" else file.path(analysis_dir, "Bonferroni_correction")

if (!file.exists(input_file)) {
  stop("Required input not found: ", input_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

loaded_objects <- load(input_file)
required_objects <- c("moduleTraitCor", "moduleTraitPvalue", "nSamples")
missing_objects <- setdiff(required_objects, loaded_objects)
if (length(missing_objects) > 0L) {
  stop("The .RData file is missing: ", paste(missing_objects, collapse = ", "))
}

if (!identical(dim(moduleTraitCor), dim(moduleTraitPvalue))) {
  stop("Correlation and P-value matrices have different dimensions.")
}
if (is.null(rownames(moduleTraitCor)) || is.null(colnames(moduleTraitCor))) {
  stop("The module-trait matrix must have row and column names.")
}

## Conservative family requested by the reviewer: 7 modules x 2 displayed traits.
n_displayed_tests <- length(moduleTraitPvalue)
if (n_displayed_tests != 14L) {
  warning("Expected 14 displayed tests, but found ", n_displayed_tests, ".")
}

moduleTraitPadjBonf <- matrix(
  p.adjust(as.vector(moduleTraitPvalue), method = "bonferroni"),
  nrow = nrow(moduleTraitPvalue),
  ncol = ncol(moduleTraitPvalue),
  dimnames = dimnames(moduleTraitPvalue)
)

moduleTraitQBH <- matrix(
  p.adjust(as.vector(moduleTraitPvalue), method = "BH"),
  nrow = nrow(moduleTraitPvalue),
  ncol = ncol(moduleTraitPvalue),
  dimnames = dimnames(moduleTraitPvalue)
)

## NC and PD are complementary indicators. This check documents the dependence
## while retaining the more conservative 14-cell correction in the main result.
complementary_traits <- FALSE
if (all(c("NC", "PD") %in% colnames(moduleTraitCor))) {
  complementary_traits <- isTRUE(all.equal(
    moduleTraitCor[, "NC"],
    -moduleTraitCor[, "PD"],
    tolerance = 1e-12
  )) && isTRUE(all.equal(
    moduleTraitPvalue[, "NC"],
    moduleTraitPvalue[, "PD"],
    tolerance = 1e-12
  ))
}

result_long <- do.call(
  rbind,
  lapply(seq_len(nrow(moduleTraitCor)), function(i) {
    data.frame(
      Module = sub("^ME", "", rownames(moduleTraitCor)[i]),
      Trait = colnames(moduleTraitCor),
      Correlation_r = as.numeric(moduleTraitCor[i, ]),
      P_nominal = as.numeric(moduleTraitPvalue[i, ]),
      P_Bonferroni_14 = as.numeric(moduleTraitPadjBonf[i, ]),
      Q_BH_14 = as.numeric(moduleTraitQBH[i, ]),
      stringsAsFactors = FALSE
    )
  })
)

result_long$Pass_Bonferroni_0.05 <- result_long$P_Bonferroni_14 < 0.05
result_long$Pass_effect_size_abs_r_gt_0.5 <- abs(result_long$Correlation_r) > 0.5
result_long$Selected_for_downstream <-
  result_long$Pass_Bonferroni_0.05 &
  result_long$Pass_effect_size_abs_r_gt_0.5

write.csv(
  result_long,
  file.path(output_dir, "WGCNA_module_trait_multiple_testing_results.csv"),
  row.names = FALSE
)

pd_results <- result_long[result_long$Trait == "PD", ]
pd_results <- pd_results[order(pd_results$P_Bonferroni_14), ]
write.csv(
  pd_results,
  file.path(output_dir, "WGCNA_PD_module_selection_results.csv"),
  row.names = FALSE
)

format_adjusted_p <- function(x) {
  ifelse(
    x < 0.001,
    formatC(x, format = "e", digits = 1),
    formatC(x, format = "f", digits = 3)
  )
}

draw_corrected_heatmap <- function() {
  module_order <- c(
    "MEturquoise", "MEbrown", "MEblue", "MEred",
    "MEgreen", "MEyellow", "MEgrey"
  )
  module_order <- module_order[module_order %in% rownames(moduleTraitCor)]
  trait_order <- c("NC", "PD")
  trait_order <- trait_order[trait_order %in% colnames(moduleTraitCor)]

  cor_mat <- moduleTraitCor[module_order, trait_order, drop = FALSE]
  padj_mat <- moduleTraitPadjBonf[module_order, trait_order, drop = FALSE]

  n_mod <- nrow(cor_mat)
  n_trait <- ncol(cor_mat)
  y_pos <- rev(seq_len(n_mod))
  palette <- colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(201)
  color_index <- function(value) {
    pmax(1L, pmin(201L, round((value + 1) / 2 * 200) + 1L))
  }

  par(mar = c(5.5, 8.0, 3.5, 4.5), family = "sans", xpd = NA)
  plot(
    NA,
    xlim = c(0.18, n_trait + 1.20),
    ylim = c(0.45, n_mod + 0.55),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Module-trait relationships\n(Bonferroni-adjusted P values)"
  )

  module_names <- sub("^ME", "", module_order)
  module_colors <- ifelse(module_names == "grey", "grey70", module_names)

  for (i in seq_len(n_mod)) {
    y <- y_pos[i]
    rect(0.20, y - 0.48, 0.36, y + 0.48,
         col = module_colors[i], border = "white")
    text(0.14, y, labels = module_order[i], adj = 1, cex = 0.92)
    for (j in seq_len(n_trait)) {
      value <- cor_mat[i, j]
      rect(j - 0.48, y - 0.48, j + 0.48, y + 0.48,
           col = palette[color_index(value)], border = "#F0F0F0")
      label <- paste0(
        formatC(value, format = "f", digits = 2),
        "\nBonf. P\n", format_adjusted_p(padj_mat[i, j])
      )
      text(j, y, labels = label, cex = 0.70,
           col = if (abs(value) >= 0.70) "white" else "black")
    }
  }

  axis(1, at = seq_len(n_trait), labels = trait_order, tick = FALSE, line = -0.5)

  ## Compact correlation color key.
  key_x0 <- n_trait + 0.68
  key_x1 <- n_trait + 0.84
  key_y <- seq(0.65, n_mod + 0.35, length.out = 201)
  for (k in 1:200) {
    rect(key_x0, key_y[k], key_x1, key_y[k + 1],
         col = palette[k], border = NA)
  }
  rect(key_x0, min(key_y), key_x1, max(key_y), border = "grey30")
  text(key_x1 + 0.05, c(min(key_y), mean(range(key_y)), max(key_y)),
       labels = c("-1", "0", "1"), adj = 0, cex = 0.75)
  text(mean(c(key_x0, key_x1)), max(key_y) + 0.20,
       labels = "r", cex = 0.85)
}

pdf(
  file.path(output_dir, "Figure3C_ModuleTrait_Bonferroni.pdf"),
  width = 5.4,
  height = 7.2,
  useDingbats = FALSE
)
draw_corrected_heatmap()
dev.off()

png(
  file.path(output_dir, "Figure3C_ModuleTrait_Bonferroni.png"),
  width = 1800,
  height = 2400,
  res = 300,
  bg = "white"
)
draw_corrected_heatmap()
dev.off()

tiff(
  file.path(output_dir, "Figure3C_ModuleTrait_Bonferroni.tif"),
  width = 1800,
  height = 2400,
  res = 300,
  compression = "lzw",
  bg = "white"
)
draw_corrected_heatmap()
dev.off()

selected_pd <- pd_results[pd_results$Selected_for_downstream, ]
summary_lines <- c(
  paste0("Input: ", if (running_from_analysis_dir) ".RData (current analysis directory)" else input_file),
  paste0("Samples used by saved WGCNA analysis: n = ", nSamples),
  paste0("Displayed tests corrected: ", n_displayed_tests),
  paste0("Bonferroni nominal-P threshold: ", format(0.05 / n_displayed_tests, digits = 6)),
  paste0("NC/PD complementary check: ", complementary_traits),
  "Selection rule: abs(r) > 0.5 and Bonferroni-adjusted P < 0.05",
  paste0(
    "Selected PD-associated modules: ",
    paste(selected_pd$Module, collapse = ", ")
  ),
  "",
  "This correction changes only multiplicity-adjusted inference; it does not alter",
  "the expression matrix, network construction, module assignments, or eigengenes."
)
writeLines(summary_lines, file.path(output_dir, "README_results.txt"))

session_file <- file.path(output_dir, "sessionInfo.txt")
sink(session_file)
cat("User library prepended: ", user_library, "\n", sep = "")
cat("External package dependencies: none (base R only)\n\n")
print(sessionInfo())
sink()

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("Outputs written to: ", output_dir, "\n", sep = "")
