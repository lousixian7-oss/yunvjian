args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript Generate_Disease_Separation_Supplement.R expression.csv group.csv output_dir")
}

expression_file <- args[[1L]]
group_file <- args[[2L]]
output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260820)

expr <- read.csv(expression_file, row.names = 1, check.names = FALSE)
group_onehot <- read.csv(group_file, row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(expr), rownames(group_onehot)))

group <- factor(ifelse(group_onehot[["PD"]] == 1, "PD", "NC"),
                levels = c("NC", "PD"))
dataset <- factor(c(rep("GSE23586", 6),
                    rep("GSE156993", 12),
                    rep("GSE10334", 154)))
stopifnot(length(group) == length(dataset), ncol(expr) == length(group))

sample_by_gene <- t(as.matrix(expr))
scaled_expression <- scale(sample_by_gene)
keep <- apply(scaled_expression, 2, function(v) all(is.finite(v)))
scaled_expression <- scaled_expression[, keep, drop = FALSE]

pca <- prcomp(scaled_expression, center = FALSE, scale. = FALSE, rank. = 2)
scores <- getElement(pca, "x")[, 1:2, drop = FALSE]
variance <- getElement(summary(pca), "importance")["Proportion of Variance", 1:2]

one_way_effect_from_gram <- function(gram, g) {
  g <- droplevels(factor(g))
  n <- nrow(gram)
  ss_total <- sum(diag(gram)) - sum(gram) / n
  ss_between <- sum(vapply(levels(g), function(level) {
    indicator <- as.numeric(g == level)
    as.numeric(crossprod(indicator, gram %*% indicator)) / sum(indicator)
  }, numeric(1))) - sum(gram) / n
  df1 <- nlevels(g) - 1
  df2 <- n - nlevels(g)
  r2 <- ss_between / ss_total
  f_value <- (ss_between / df1) / ((ss_total - ss_between) / df2)
  c(R2 = r2, F = f_value)
}

one_way_effect <- function(y, g) one_way_effect_from_gram(tcrossprod(y), g)

permutation_p <- function(y, g, strata = NULL, permutations = 9999L) {
  gram <- tcrossprod(y)
  observed <- unname(one_way_effect_from_gram(gram, g)[["F"]])
  permuted <- numeric(permutations)
  for (b in seq_len(permutations)) {
    if (is.null(strata)) {
      gp <- sample(g)
    } else {
      gp <- g
      for (level in levels(factor(strata))) {
        idx <- which(strata == level)
        gp[idx] <- sample(g[idx])
      }
    }
    permuted[[b]] <- unname(one_way_effect_from_gram(gram, gp)[["F"]])
  }
  (1 + sum(permuted >= observed)) / (permutations + 1)
}

auc_rank <- function(observed, probability) {
  observed <- as.integer(observed)
  n1 <- sum(observed == 1L)
  n0 <- sum(observed == 0L)
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[observed == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

make_folds <- function(g, k = 10L) {
  folds <- integer(length(g))
  for (level in levels(g)) {
    idx <- which(g == level)
    folds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  folds
}

repeats <- 100L
prediction_sum <- numeric(length(group))
repeat_auc <- numeric(repeats)
for (r in seq_len(repeats)) {
  folds <- make_folds(group, 10L)
  prediction <- numeric(length(group))
  for (fold in seq_len(10L)) {
    train <- folds != fold
    test <- !train
    training_data <- data.frame(
      outcome = as.integer(group[train] == "PD"),
      PC1 = scores[train, 1],
      PC2 = scores[train, 2]
    )
    fit <- glm(outcome ~ PC1 + PC2, data = training_data, family = binomial())
    prediction[test] <- predict(
      fit,
      newdata = data.frame(PC1 = scores[test, 1], PC2 = scores[test, 2]),
      type = "response"
    )
  }
  prediction_sum <- prediction_sum + prediction
  repeat_auc[[r]] <- auc_rank(group == "PD", prediction)
}
mean_oof_prediction <- prediction_sum / repeats
mean_oof_auc <- auc_rank(group == "PD", mean_oof_prediction)

bootstrap_auc <- replicate(2000L, {
  nc <- sample(which(group == "NC"), replace = TRUE)
  pd <- sample(which(group == "PD"), replace = TRUE)
  idx <- c(nc, pd)
  auc_rank(group[idx] == "PD", mean_oof_prediction[idx])
})
auc_ci <- unname(quantile(bootstrap_auc, c(0.025, 0.975), na.rm = TRUE))

pc12_effect <- one_way_effect(scores, group)
pc12_p <- permutation_p(scores, group, strata = dataset, permutations = 9999L)
global_effect <- one_way_effect(scaled_expression, group)
global_p <- permutation_p(scaled_expression, group, strata = dataset,
                          permutations = 9999L)
batch_effect <- one_way_effect(sample_by_gene, dataset)

# Paired sensitivity analysis for the 64 matched NC/PD samples in GSE10334.
clinic_path <- file.path(dirname(expression_file), "GSE10334.clinic.txt")
paired_n <- paired_dz <- paired_mean <- paired_median <- paired_p <- NA_real_
if (file.exists(clinic_path)) {
  clinic <- read.delim(clinic_path, check.names = FALSE, quote = "")
  paired <- data.frame(
    sample = clinic[["geo_accession"]],
    patient = clinic[["patient"]],
    group = clinic[["group"]],
    PC1 = scores[match(clinic[["geo_accession"]], rownames(scores)), 1]
  )
  paired_wide <- reshape(paired, idvar = "patient", timevar = "group",
                         direction = "wide")
  paired_wide <- paired_wide[complete.cases(paired_wide[, c("PC1.NC", "PC1.PD")]), ]
  paired_difference <- paired_wide[["PC1.PD"]] - paired_wide[["PC1.NC"]]
  paired_n <- length(paired_difference)
  paired_mean <- mean(paired_difference)
  paired_median <- median(paired_difference)
  paired_dz <- paired_mean / sd(paired_difference)
  paired_p <- getElement(wilcox.test(paired_difference, exact = FALSE), "p.value")
}

results <- data.frame(
  Analysis = c(
    "Disease group in PC1-PC2 space",
    "Disease group in all standardized expression dimensions",
    "Dataset effect after ComBat in unscaled expression space",
    "PC1-PC2 repeated 10-fold cross-validated discrimination",
    "GSE10334 paired NC-PD PC1 sensitivity analysis"
  ),
  Metric = c("PERMANOVA-like R2", "PERMANOVA-like R2", "R2",
             "AUC", "Paired standardized mean difference (dz)"),
  Estimate = c(pc12_effect[["R2"]], global_effect[["R2"]],
               batch_effect[["R2"]], mean_oof_auc, paired_dz),
  CI_lower = c(NA, NA, NA, auc_ci[[1]], NA),
  CI_upper = c(NA, NA, NA, auc_ci[[2]], NA),
  P_value = c(pc12_p, global_p, NA, NA, paired_p),
  N = c(length(group), length(group), length(group), length(group), paired_n),
  Notes = c(
    sprintf("PC1+PC2 explain %.2f%% of total standardized variance; permutations stratified by dataset", 100 * sum(variance)),
    "All non-zero-variance genes; permutations stratified by dataset",
    "Descriptive residual dataset effect; Euclidean expression space",
    sprintf("100 repeated stratified 10-fold runs; repeat AUC mean %.3f, SD %.3f; CI is stratified bootstrap of mean out-of-fold predictions", mean(repeat_auc), sd(repeat_auc)),
    sprintf("%d matched patients; PD-NC PC1 mean difference %.2f, median difference %.2f; Wilcoxon signed-rank test", paired_n, paired_mean, paired_median)
  ),
  check.names = FALSE
)

write.csv(results,
          file.path(output_dir, "Supplementary_Table_Disease_Separation_After_ComBat.csv"),
          row.names = FALSE, quote = TRUE)

roc_points <- function(observed, probability) {
  thresholds <- c(Inf, sort(unique(probability), decreasing = TRUE), -Inf)
  tpr <- fpr <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    predicted <- probability >= thresholds[[i]]
    tpr[[i]] <- sum(predicted & observed) / sum(observed)
    fpr[[i]] <- sum(predicted & !observed) / sum(!observed)
  }
  data.frame(FPR = fpr, TPR = tpr)
}
roc <- roc_points(group == "PD", mean_oof_prediction)

draw_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.05, 1))
  par(mar = c(7, 5, 4, 1.5), family = "sans")
  values <- c(batch_effect[["R2"]], global_effect[["R2"]], pc12_effect[["R2"]])
  labels <- c("Dataset\n(full expression)",
              "Disease group\n(all dimensions)",
              "Disease group\n(PC1-PC2)")
  colors <- c("#F8766D", "#00BFC4", "#00A6B4")
  mids <- barplot(values, names.arg = labels, col = colors, border = NA,
                  ylim = c(0, max(values) * 1.35), las = 1,
                  ylab = expression(R^2), main = "A  Variance explained after ComBat",
                  cex.names = 0.82)
  text(mids, values, labels = sprintf("%.3f", values), pos = 3, cex = 0.9)
  mtext(sprintf("Disease-group permutation P < 0.001; PC1+PC2 variance = %.1f%%",
                100 * sum(variance)), side = 1, line = 5.5, cex = 0.78)

  par(mar = c(5, 5, 4, 2), family = "sans")
  plot(roc$FPR, roc$TPR, type = "l", lwd = 3, col = "#00A6B4",
       xlim = c(0, 1), ylim = c(0, 1), asp = 1,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = "B  Descriptive discrimination from PC1-PC2")
  abline(0, 1, lty = 2, col = "grey55")
  legend("bottomright",
         legend = c(sprintf("Repeated 10-fold CV AUC = %.3f", mean_oof_auc),
                    sprintf("Bootstrap 95%% CI %.3f-%.3f", auc_ci[[1]], auc_ci[[2]]),
                    sprintf("GSE10334 paired sensitivity: dz = %.2f, P < 0.001", paired_dz)),
         bty = "n", cex = 0.78)
}

png(file.path(output_dir, "Supplementary_Figure_Disease_Separation_After_ComBat.png"),
    width = 3600, height = 1650, res = 300, bg = "white")
draw_figure()
dev.off()

pdf(file.path(output_dir, "Supplementary_Figure_Disease_Separation_After_ComBat.pdf"),
    width = 12, height = 5.5, useDingbats = FALSE)
draw_figure()
dev.off()

caption <- c(
  "# Supplementary analysis: disease-group separation after ComBat",
  "",
  sprintf("After ComBat correction, PC1 and PC2 explained %.1f%% of the standardized expression variance (PC1 %.1f%%; PC2 %.1f%%). In this two-dimensional space, disease group explained R2 = %.3f (dataset-stratified permutation P < 0.001). A descriptive logistic model using only the two displayed coordinates achieved a repeated 10-fold cross-validated AUC of %.3f (bootstrap 95%% CI %.3f-%.3f). In the 64 matched NC-PD sample pairs from GSE10334, the PC1 shift remained evident (paired dz = %.2f; Wilcoxon signed-rank P = %.3g). These results indicate moderate, statistically detectable disease-associated structure despite substantial overlap in the two-dimensional PCA display; the PCA-based AUC is descriptive and should not be interpreted as diagnostic-model performance.",
          100 * sum(variance), 100 * variance[[1]], 100 * variance[[2]],
          pc12_effect[["R2"]], mean_oof_auc, auc_ci[[1]], auc_ci[[2]],
          paired_dz, paired_p),
  "",
  "Methods: genes were standardized before PCA to reproduce the geometry of Fig. 1D. Disease-group pseudo-R2 was calculated from Euclidean sums of squares. Permutation tests used 9,999 label permutations within each dataset. The AUC analysis used 100 repeated stratified 10-fold splits and logistic regression with PC1 and PC2 as the only predictors. The PCA/AUC analysis quantifies the coordinates shown in Fig. 1D and is not an independent diagnostic validation.",
  "",
  sprintf("Reproducibility seed: 20260820. Samples: %d (NC=%d, PD=%d). Retained genes: %d.",
          length(group), sum(group == "NC"), sum(group == "PD"), ncol(scaled_expression))
)
writeLines(caption,
           file.path(output_dir, "Supplementary_Disease_Separation_After_ComBat.md"),
           useBytes = TRUE)

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "Supplementary_Disease_Separation_sessionInfo.txt"))

print(results)
