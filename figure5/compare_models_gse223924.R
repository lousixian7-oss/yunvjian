options(stringsAsFactors = FALSE)
set.seed(20260902)

user_library <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
.libPaths(c(user_library, .libPaths()))

input_dir <- "model_revision_inputs"
gse_dir <- "GSE223924_cross_platform_validation"
out_dir <- "three_vs_alt4_model_sensitivity_20260902"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_sets <- list(
  Locked4_C3 = c("AGT", "C3", "CXCR4", "FOS"),
  Three_gene = c("AGT", "CXCR4", "FOS"),
  Alt4_EPHX2 = c("AGT", "EPHX2", "CXCR4", "FOS")
)
all_genes <- unique(unlist(model_sets))

auc_rank <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
boot_auc <- function(y, p, B = 5000, seed = 1) {
  set.seed(seed); i0 <- which(y == 0); i1 <- which(y == 1); z <- numeric(B)
  for (b in seq_len(B)) {
    ii <- c(sample(i0, length(i0), replace = TRUE), sample(i1, length(i1), replace = TRUE))
    z[b] <- auc_rank(y[ii], p[ii])
  }
  unname(quantile(z, c(.025, .975), names = FALSE))
}
paired_delta <- function(y, p_new, p_ref, B = 5000, seed = 1) {
  set.seed(seed); i0 <- which(y == 0); i1 <- which(y == 1); z <- numeric(B)
  for (b in seq_len(B)) {
    ii <- c(sample(i0, length(i0), replace = TRUE), sample(i1, length(i1), replace = TRUE))
    z[b] <- auc_rank(y[ii], p_new[ii]) - auc_rank(y[ii], p_ref[ii])
  }
  c(delta = auc_rank(y, p_new) - auc_rank(y, p_ref),
    low = quantile(z, .025, names = FALSE), high = quantile(z, .975, names = FALSE),
    p = min(1, 2 * min(mean(z <= 0), mean(z >= 0))))
}

message("[1/5] Reading GSE223924 and fixed development training split...")
raw <- read.csv(gzfile(file.path(gse_dir, "GSE223924_raw_count.csv.gz")), row.names = 1, check.names = FALSE)
lcpm <- read.csv(gzfile(file.path(gse_dir, "GSE223924_lcpm_exp.csv.gz")), row.names = 1, check.names = FALSE)
dev_expr <- read.csv(file.path(input_dir, "development_expression.csv"), row.names = 1, check.names = FALSE)
dev_info <- read.csv(file.path(input_dir, "development_sample_info.csv"), row.names = 1, check.names = FALSE)
split <- read.csv(file.path(input_dir, "patient_split.csv"), check.names = FALSE)
stopifnot(all(all_genes %in% rownames(raw)), all(all_genes %in% rownames(lcpm)), all(all_genes %in% rownames(dev_expr)))

target_ids <- c(grep("^H[0-9]+$", colnames(raw), value = TRUE), grep("^PT[0-9]+$", colnames(raw), value = TRUE))
stopifnot(length(target_ids) == 20L, sum(grepl("^H", target_ids)) == 10L, sum(grepl("^PT", target_ids)) == 10L)
meta <- data.frame(sample_id = target_ids, group = ifelse(grepl("^PT", target_ids), "PD", "NC"), stringsAsFactors = FALSE)
meta$y <- as.integer(meta$group == "PD")

positive <- rowSums(raw > 0) == ncol(raw)
geo <- exp(rowMeans(log(as.matrix(raw[positive, , drop = FALSE]))))
ratios <- sweep(as.matrix(raw[positive, , drop = FALSE]), 1, geo, "/")
sf <- apply(ratios, 2, median, na.rm = TRUE)
sf <- sf / exp(mean(log(sf)))
log_raw <- log2(sweep(as.matrix(raw), 2, sf, "/") + 1)
raw_target <- t(log_raw[all_genes, target_ids, drop = FALSE])
lcpm_target <- t(as.matrix(lcpm[all_genes, target_ids, drop = FALSE]))

train_ids <- split$sample_id[split$cohort == "Training"]
train_all <- as.data.frame(t(dev_expr[all_genes, train_ids, drop = FALSE]), check.names = FALSE)
train_y <- as.integer(dev_info[train_ids, "title"] == "Treat")

message("[2/5] Fitting training-only models and producing label-free target-cohort predictions...")
fits <- list(); pred_raw <- list(); pred_lcpm <- list()
for (nm in names(model_sets)) {
  genes <- model_sets[[nm]]
  mu <- vapply(train_all[, genes, drop = FALSE], mean, numeric(1))
  sig <- vapply(train_all[, genes, drop = FALSE], sd, numeric(1))
  train_z <- as.data.frame(sweep(sweep(as.matrix(train_all[, genes, drop = FALSE]), 2, mu, "-"), 2, sig, "/"), check.names = FALSE)
  fits[[nm]] <- glm(train_y ~ ., data = train_z, family = binomial())
  raw_z <- as.data.frame(scale(raw_target[, genes, drop = FALSE]), check.names = FALSE)
  lcpm_z <- as.data.frame(scale(lcpm_target[, genes, drop = FALSE]), check.names = FALSE)
  pred_raw[[nm]] <- predict(fits[[nm]], newdata = raw_z, type = "response")
  pred_lcpm[[nm]] <- predict(fits[[nm]], newdata = lcpm_z, type = "response")
}

message("[3/5] Computing AUCs and bootstrap intervals...")
metric_rows <- list(); q <- 1L
for (prep in c("Median_of_ratios_log2", "Official_log2_CPM")) for (nm in names(model_sets)) {
  p <- if (prep == "Median_of_ratios_log2") pred_raw[[nm]] else pred_lcpm[[nm]]
  ci <- boot_auc(meta$y, p, seed = 6100 + q)
  metric_rows[[q]] <- data.frame(
    model = nm, genes = paste(model_sets[[nm]], collapse = ";"), preprocessing = prep,
    samples = 20L, NC = 10L, PD = 10L, AUC = auc_rank(meta$y, p),
    CI_low = ci[1], CI_high = ci[2], Brier = mean((meta$y - p)^2)
  )
  q <- q + 1L
}
metrics <- do.call(rbind, metric_rows)

delta_rows <- list(); q <- 1L
for (prep in c("Median_of_ratios_log2", "Official_log2_CPM")) for (nm in c("Three_gene", "Alt4_EPHX2")) {
  preds <- if (prep == "Median_of_ratios_log2") pred_raw else pred_lcpm
  z <- paired_delta(meta$y, preds[[nm]], preds$Locked4_C3, seed = 7100 + q)
  delta_rows[[q]] <- data.frame(
    comparison = paste0(nm, " minus Locked4_C3"), preprocessing = prep,
    delta_AUC = z["delta"], CI_low = z["low"], CI_high = z["high"],
    bootstrap_two_sided_P = z["p"]
  )
  q <- q + 1L
}
deltas <- do.call(rbind, delta_rows)

message("[4/5] Writing outputs...")
predictions <- meta
for (nm in names(model_sets)) {
  predictions[[paste0(nm, "_probability_MOR")]] <- pred_raw[[nm]]
  predictions[[paste0(nm, "_probability_lcpm")]] <- pred_lcpm[[nm]]
}
write.csv(metrics, file.path(out_dir, "Table_GSE223924_three_model_AUC.csv"), row.names = FALSE)
write.csv(deltas, file.path(out_dir, "Table_GSE223924_paired_AUC_differences.csv"), row.names = FALSE)
write.csv(predictions, file.path(out_dir, "GSE223924_three_model_predictions.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo_GSE223924.txt"))

message("[5/5] Complete.")
print(metrics)
print(deltas)
