options(stringsAsFactors = FALSE)
set.seed(20260902)

user_library <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
.libPaths(c(user_library, .libPaths()))
required_packages <- c("pROC")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))

input_dir <- "model_revision_inputs"
out_dir <- "three_vs_alt4_model_sensitivity_20260902"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_sets <- list(
  Locked4_C3 = c("AGT", "C3", "CXCR4", "FOS"),
  Three_gene = c("AGT", "CXCR4", "FOS"),
  Alt4_EPHX2 = c("AGT", "EPHX2", "CXCR4", "FOS")
)
all_genes <- unique(unlist(model_sets))

message("[1/7] Reading fixed development split and external cohorts...")
dev_expr <- read.csv(file.path(input_dir, "development_expression.csv"), row.names = 1, check.names = FALSE)
dev_info <- read.csv(file.path(input_dir, "development_sample_info.csv"), row.names = 1, check.names = FALSE)
split_info <- read.csv(file.path(input_dir, "patient_split.csv"), check.names = FALSE)
same_expr <- read.csv(file.path(input_dir, "GSE16134_expression.csv"), row.names = 1, check.names = FALSE)
same_info <- read.csv(file.path(input_dir, "GSE16134_sample_info.csv"), row.names = 1, check.names = FALSE)
counts <- read.delim(gzfile(file.path(input_dir, "GSE173078_counts.txt.gz")), check.names = FALSE)
if (!all(all_genes %in% rownames(dev_expr)) || !all(all_genes %in% rownames(same_expr)) ||
    !all(all_genes %in% rownames(counts))) stop("At least one candidate gene is absent from an analysis dataset.")

dev <- as.data.frame(t(dev_expr[all_genes, , drop = FALSE]), check.names = FALSE)
dev$sample_id <- rownames(dev)
dev <- merge(dev, split_info[, c("sample_id", "patient_id", "dataset", "cohort")],
             by = "sample_id", sort = FALSE)
dev$y <- as.integer(dev_info[dev$sample_id, "title"] == "Treat")
train <- dev[dev$cohort == "Training", , drop = FALSE]
test <- dev[dev$cohort == "Internal test", , drop = FALSE]
stopifnot(nrow(train) == 120L, nrow(test) == 52L)
stopifnot(length(intersect(unique(train$patient_id), unique(test$patient_id))) == 0L)

same <- as.data.frame(t(same_expr[all_genes, , drop = FALSE]), check.names = FALSE)
same$sample_id <- rownames(same)
same$y <- as.integer(same_info[same$sample_id, "Group"] == "PD")
same$patient_num <- sub(".*[Pp]atient[ ]*([0-9]+).*", "\\1", same_info[same$sample_id, "Title"])
same$patient_id <- paste0("GSE16134_P", same$patient_num)
same$dataset <- "GSE16134"

positive <- rowSums(counts > 0) == ncol(counts)
geo_mean <- exp(rowMeans(log(as.matrix(counts[positive, , drop = FALSE]))))
ratios <- sweep(as.matrix(counts[positive, , drop = FALSE]), 1, geo_mean, "/")
size_factor <- apply(ratios, 2, median, na.rm = TRUE)
size_factor <- size_factor / exp(mean(log(size_factor)))
log_expr <- log2(sweep(as.matrix(counts), 2, size_factor, "/") + 1)
series_lines <- readLines(gzfile(file.path(input_dir, "GSE173078_series_matrix.txt.gz")), warn = FALSE)
title_line <- grep("^!Sample_title", series_lines, value = TRUE)[1]
titles <- gsub('^"|"$', "", strsplit(title_line, "\t", fixed = TRUE)[[1]][-1])
rna_meta <- data.frame(sample_id = sub("_.*", "", titles), title = titles, stringsAsFactors = FALSE)
rna_meta$phenotype <- ifelse(grepl("Periodontitis", rna_meta$title, ignore.case = TRUE), "PD",
                             ifelse(grepl("Healthy", rna_meta$title, ignore.case = TRUE), "NC", "Gingivitis"))
rna <- rna_meta[rna_meta$phenotype %in% c("NC", "PD"), , drop = FALSE]
rna$y <- as.integer(rna$phenotype == "PD")
rna$patient_id <- paste0("GSE173078_", rna$sample_id)
rna$dataset <- "GSE173078"
rna <- cbind(rna, as.data.frame(t(log_expr[all_genes, rna$sample_id, drop = FALSE]), check.names = FALSE))

auc_rank <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
z_apply <- function(x, mu, sig) sweep(sweep(as.matrix(x), 2, mu, "-"), 2, sig, "/")
cohort_z <- function(x) scale(as.matrix(x))

make_group_folds <- function(d, k = 5, seed = 1) {
  set.seed(seed)
  ids <- unique(d$patient_id)
  patient_frame <- do.call(rbind, lapply(ids, function(id) {
    z <- d[d$patient_id == id, , drop = FALSE]
    signature <- if (length(unique(z$y)) == 2) "Mixed" else if (z$y[1] == 1) "PD" else "NC"
    data.frame(patient_id = id, stratum = paste(z$dataset[1], signature, sep = "__"))
  }))
  patient_frame$fold <- NA_integer_
  for (st in unique(patient_frame$stratum)) {
    ii <- which(patient_frame$stratum == st)
    patient_frame$fold[ii] <- sample(rep(seq_len(k), length.out = length(ii)))
  }
  setNames(patient_frame$fold, patient_frame$patient_id)
}

message("[2/7] Running 20x5-fold patient-grouped training CV...")
cv_rows <- list(); rr <- 1L
for (replicate_id in seq_len(20)) {
  folds <- make_group_folds(train, k = 5, seed = 1000 + replicate_id)
  for (fold_id in seq_len(5)) {
    va <- folds[train$patient_id] == fold_id
    tr <- !va
    for (model_name in names(model_sets)) {
      genes <- model_sets[[model_name]]
      mu <- vapply(train[tr, genes, drop = FALSE], mean, numeric(1))
      sig <- vapply(train[tr, genes, drop = FALSE], sd, numeric(1))
      xtr <- as.data.frame(z_apply(train[tr, genes, drop = FALSE], mu, sig), check.names = FALSE)
      xva <- as.data.frame(z_apply(train[va, genes, drop = FALSE], mu, sig), check.names = FALSE)
      fit <- suppressWarnings(glm(train$y[tr] ~ ., data = xtr, family = binomial()))
      pred <- suppressWarnings(predict(fit, newdata = xva, type = "response"))
      cv_rows[[rr]] <- data.frame(replicate = replicate_id, fold = fold_id,
                                  model = model_name, n = sum(va), auc = auc_rank(train$y[va], pred))
      rr <- rr + 1L
    }
  }
}
cv <- do.call(rbind, cv_rows)
cv_summary <- do.call(rbind, lapply(split(cv, cv$model), function(z) data.frame(
  model = z$model[1], genes = paste(model_sets[[z$model[1]]], collapse = ";"),
  valid_folds = sum(is.finite(z$auc)), mean_fold_auc = mean(z$auc, na.rm = TRUE),
  sd_fold_auc = sd(z$auc, na.rm = TRUE), median_fold_auc = median(z$auc, na.rm = TRUE)
)))

message("[3/7] Fitting locked models and evaluating cohorts...")
evaluate_model <- function(model_name, genes) {
  mu <- vapply(train[, genes, drop = FALSE], mean, numeric(1))
  sig <- vapply(train[, genes, drop = FALSE], sd, numeric(1))
  tr_z <- as.data.frame(z_apply(train[, genes, drop = FALSE], mu, sig), check.names = FALSE)
  te_z <- as.data.frame(z_apply(test[, genes, drop = FALSE], mu, sig), check.names = FALSE)
  same_z <- as.data.frame(cohort_z(same[, genes, drop = FALSE]), check.names = FALSE)
  rna_z <- as.data.frame(cohort_z(rna[, genes, drop = FALSE]), check.names = FALSE)
  fit <- suppressWarnings(glm(train$y ~ ., data = tr_z, family = binomial()))
  list(name = model_name, genes = genes, fit = fit,
       train = predict(fit, newdata = tr_z, type = "response"),
       test = predict(fit, newdata = te_z, type = "response"),
       same = predict(fit, newdata = same_z, type = "response"),
       rna = predict(fit, newdata = rna_z, type = "response"))
}
fits <- lapply(names(model_sets), function(nm) evaluate_model(nm, model_sets[[nm]]))
names(fits) <- names(model_sets)
cohort_objects <- list(Training = train, Internal_test = test, GSE16134 = same, GSE173078 = rna)

metric_rows <- list(); q <- 1L
for (nm in names(fits)) for (co in names(cohort_objects)) {
  d <- cohort_objects[[co]]
  key <- switch(co, Training = "train", Internal_test = "test", GSE16134 = "same", GSE173078 = "rna")
  p <- fits[[nm]][[key]]
  metric_rows[[q]] <- data.frame(
    model = nm, genes = paste(model_sets[[nm]], collapse = ";"), cohort = co,
    samples = nrow(d), patients = length(unique(d$patient_id)),
    AUC = auc_rank(d$y, p), Brier = mean((d$y - p)^2)
  )
  q <- q + 1L
}
metrics <- do.call(rbind, metric_rows)

message("[4/7] Computing AIC, nested C3 deletion test, and coefficients...")
fit_stats <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]$fit
  data.frame(model = nm, genes = paste(model_sets[[nm]], collapse = ";"),
             logLik = as.numeric(logLik(f)), AIC = AIC(f),
             BIC_patient_penalty = -2 * as.numeric(logLik(f)) + log(length(unique(train$patient_id))) * length(coef(f)))
}))
lrt_c3 <- anova(fits$Three_gene$fit, fits$Locked4_C3$fit, test = "LRT")
lrt_table <- data.frame(
  comparison = "Three_gene versus Locked4_C3 (incremental C3)",
  deviance_difference = lrt_c3$Deviance[2], df_difference = lrt_c3$Df[2],
  likelihood_ratio_P = lrt_c3$`Pr(>Chi)`[2]
)
lrt_ephx2 <- anova(fits$Three_gene$fit, fits$Alt4_EPHX2$fit, test = "LRT")
lrt_ephx2_table <- data.frame(
  comparison = "Three_gene versus Alt4_EPHX2 (incremental EPHX2)",
  deviance_difference = lrt_ephx2$Deviance[2], df_difference = lrt_ephx2$Df[2],
  likelihood_ratio_P = lrt_ephx2$`Pr(>Chi)`[2]
)
coefficient_table <- do.call(rbind, lapply(names(fits), function(nm) {
  z <- as.data.frame(summary(fits[[nm]]$fit)$coefficients)
  z$term <- rownames(z); rownames(z) <- NULL
  names(z)[1:4] <- c("estimate", "standard_error", "z_value", "Wald_P")
  cbind(model = nm, z[, c("term", "estimate", "standard_error", "z_value", "Wald_P")])
}))

cluster_boot_auc <- function(d, p, B = 5000, seed = 1) {
  set.seed(seed); ids <- unique(d$patient_id); ans <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    picked <- sample(ids, length(ids), replace = TRUE)
    idx <- unlist(lapply(picked, function(id) which(d$patient_id == id)), use.names = FALSE)
    ans[b] <- auc_rank(d$y[idx], p[idx])
  }
  ans <- ans[is.finite(ans)]
  unname(quantile(ans, c(.025, .975), names = FALSE))
}
paired_boot_delta <- function(d, p_new, p_ref, B = 5000, seed = 1) {
  set.seed(seed); ids <- unique(d$patient_id); ans <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    picked <- sample(ids, length(ids), replace = TRUE)
    idx <- unlist(lapply(picked, function(id) which(d$patient_id == id)), use.names = FALSE)
    ans[b] <- auc_rank(d$y[idx], p_new[idx]) - auc_rank(d$y[idx], p_ref[idx])
  }
  ans <- ans[is.finite(ans)]
  c(delta = auc_rank(d$y, p_new) - auc_rank(d$y, p_ref),
    low = quantile(ans, .025, names = FALSE), high = quantile(ans, .975, names = FALSE),
    p = min(1, 2 * min(mean(ans <= 0), mean(ans >= 0))))
}

message("[5/7] Running patient-cluster bootstrap inference...")
ci_rows <- list(); q <- 1L
for (nm in names(fits)) for (co in names(cohort_objects)) {
  d <- cohort_objects[[co]]
  key <- switch(co, Training = "train", Internal_test = "test", GSE16134 = "same", GSE173078 = "rna")
  ci <- cluster_boot_auc(d, fits[[nm]][[key]], seed = 3000 + q)
  ci_rows[[q]] <- data.frame(model = nm, cohort = co, AUC = auc_rank(d$y, fits[[nm]][[key]]),
                             CI_low = ci[1], CI_high = ci[2])
  q <- q + 1L
}
ci_table <- do.call(rbind, ci_rows)

delta_rows <- list(); q <- 1L
for (nm in c("Three_gene", "Alt4_EPHX2")) for (co in c("Internal_test", "GSE16134", "GSE173078")) {
  d <- cohort_objects[[co]]
  key <- switch(co, Internal_test = "test", GSE16134 = "same", GSE173078 = "rna")
  z <- paired_boot_delta(d, fits[[nm]][[key]], fits$Locked4_C3[[key]], seed = 5000 + q)
  delta_rows[[q]] <- data.frame(
    comparison = paste0(nm, " minus Locked4_C3"), cohort = co,
    delta_AUC = z["delta"], CI_low = z["low"], CI_high = z["high"],
    bootstrap_two_sided_P = z["p"]
  )
  q <- q + 1L
}
delta_table <- do.call(rbind, delta_rows)

message("[6/7] Writing auditable outputs...")
write.csv(cv, file.path(out_dir, "Table_patient_grouped_CV_all_folds.csv"), row.names = FALSE)
write.csv(cv_summary[order(-cv_summary$mean_fold_auc), ], file.path(out_dir, "Table_patient_grouped_CV_summary.csv"), row.names = FALSE)
write.csv(metrics, file.path(out_dir, "Table_all_cohort_metrics.csv"), row.names = FALSE)
write.csv(ci_table, file.path(out_dir, "Table_AUC_cluster_bootstrap_CI.csv"), row.names = FALSE)
write.csv(delta_table, file.path(out_dir, "Table_paired_AUC_differences_vs_locked4.csv"), row.names = FALSE)
write.csv(fit_stats[order(fit_stats$AIC), ], file.path(out_dir, "Table_training_information_criteria.csv"), row.names = FALSE)
write.csv(lrt_table, file.path(out_dir, "Table_incremental_C3_LRT.csv"), row.names = FALSE)
write.csv(lrt_ephx2_table, file.path(out_dir, "Table_incremental_EPHX2_LRT.csv"), row.names = FALSE)
write.csv(coefficient_table, file.path(out_dir, "Table_standardized_model_coefficients.csv"), row.names = FALSE)
saveRDS(fits$Three_gene$fit, file.path(out_dir, "Training_model_three_gene.rds"))
saveRDS(fits$Alt4_EPHX2$fit, file.path(out_dir, "Training_model_alt4_EPHX2.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))

report <- c(
  "# Three-gene and EPHX2 four-gene sensitivity analysis", "",
  "## Design",
  "- Reference: locked AGT + C3 + CXCR4 + FOS model.",
  "- Candidate 1: AGT + CXCR4 + FOS.",
  "- Candidate 2: AGT + EPHX2 + CXCR4 + FOS.",
  "- Model fitting was restricted to the original patient-separated training cohort.",
  "- Cross-validation used 20 repeats of five-fold patient-grouped splitting (100 folds).",
  "- External cohorts were used only for locked evaluation and not model selection.", "",
  "## Key outputs",
  paste(capture.output(print(cv_summary[order(-cv_summary$mean_fold_auc), ])), collapse = "\n"), "",
  paste(capture.output(print(metrics)), collapse = "\n"), "",
  paste(capture.output(print(delta_table)), collapse = "\n"), "",
  paste0("- Incremental C3 likelihood-ratio P = ", signif(lrt_table$likelihood_ratio_P, 4)), "",
  paste0("- Incremental EPHX2 likelihood-ratio P = ", signif(lrt_ephx2_table$likelihood_ratio_P, 4)), "",
  "## Interpretation rule",
  "- Prefer parsimony only if the simpler model does not materially reduce patient-grouped CV or locked validation performance.",
  "- Do not select among these models using an external-cohort point estimate alone.",
  "- GSE173078 is a small cross-platform cohort and remains exploratory."
)
writeLines(report, file.path(out_dir, "Comparison_report.md"), useBytes = TRUE)

message("[7/7] Complete.")
print(cv_summary[order(-cv_summary$mean_fold_auc), ])
print(metrics)
print(fit_stats[order(fit_stats$AIC), ])
print(lrt_table)
print(lrt_ephx2_table)
print(delta_table)
