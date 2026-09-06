options(stringsAsFactors = FALSE)
set.seed(20260902)

## Reuse the user's existing R libraries; do not install packages here.
.libPaths(c(
  "C:/Users/32266/AppData/Local/R/win-library/4.5",
  "G:/gurobi/Rlib45",
  .libPaths()
))

required_packages <- c("ggplot2", "magick")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing existing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2)
  library(magick)
})

base_dir <- "."
input_dir <- file.path(base_dir, "model_revision_inputs")
sensitivity_dir <- file.path(base_dir, "three_vs_alt4_model_sensitivity_20260902")
out_dir <- file.path(base_dir, "three_gene_Figure5_FigureS2_20260902")
figure_dir <- out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "panels"), recursive = TRUE, showWarnings = FALSE)

genes <- c("AGT", "CXCR4", "FOS")
group_cols <- c(NC = "#27C4CC", PD = "#F47B72")
model_col <- "#EE6A5F"
cohort_cols <- c(Training = "#D95F5F", `Internal test` = "#2F7FBF", GSE16134 = "#22A884", GSE223924 = "#E08214")
gene_cols <- c(AGT = "#E6B566", CXCR4 = "#8B80D7", FOS = "#58A6D8", `Three-gene model` = "#D73027")

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.05)),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = rel(.82)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey55"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.spacing = grid::unit(.75, "lines"),
      plot.margin = margin(9, 10, 8, 10)
    )
}

panel_png <- function(name, p, width, height, dpi = 210) {
  f <- file.path(out_dir, "panels", paste0(name, ".png"))
  ggsave(f, p, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(out_dir, "panels", paste0(name, ".pdf")), p,
         width = width, height = height, device = cairo_pdf, bg = "white")
  f
}

auc_rank <- function(y, p) {
  y <- as.integer(y); p <- as.numeric(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_df <- function(y, p, cohort) {
  ord <- order(p, decreasing = TRUE); yy <- y[ord]
  data.frame(FPR = c(0, cumsum(yy == 0) / sum(yy == 0), 1),
             TPR = c(0, cumsum(yy == 1) / sum(yy == 1), 1), Cohort = cohort)
}

calibration_stats <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6); lp <- qlogis(p)
  slope_fit <- suppressWarnings(glm(y ~ lp, family = binomial()))
  int_fit <- suppressWarnings(glm(y ~ 1, offset = lp, family = binomial()))
  c(Brier = mean((y - p)^2), Intercept = unname(coef(int_fit)[1]), Slope = unname(coef(slope_fit)[2]))
}

calibration_bins <- function(d, cohort, bins = 6) {
  q <- unique(quantile(d$probability, seq(0, 1, length.out = bins + 1), na.rm = TRUE))
  d$bin <- cut(d$probability, breaks = q, include.lowest = TRUE, labels = FALSE)
  z <- aggregate(cbind(probability, y) ~ bin, d, mean); z$Cohort <- cohort; z
}

net_benefit <- function(y, p, thresholds) {
  prev <- mean(y == 1)
  model <- vapply(thresholds, function(t) mean(p >= t & y == 1) - mean(p >= t & y == 0) * t / (1 - t), numeric(1)) / prev
  all <- (prev - (1 - prev) * thresholds / (1 - thresholds)) / prev
  data.frame(Threshold = thresholds, Model = model, All = all, None = 0)
}

patient_test <- function(d, value_col) {
  z <- aggregate(d[[value_col]], list(patient_id = d$patient_id, group = d$group), mean)
  names(z)[3] <- "value"
  w <- reshape(z, idvar = "patient_id", timevar = "group", direction = "wide")
  if (all(c("value.NC", "value.PD") %in% names(w))) {
    paired <- w[complete.cases(w[, c("value.NC", "value.PD")]), , drop = FALSE]
  } else paired <- w[0, , drop = FALSE]
  if (nrow(paired) >= 5) {
    p <- suppressWarnings(wilcox.test(paired$value.PD, paired$value.NC, paired = TRUE, exact = FALSE)$p.value)
    method <- paste0("paired patient-level Wilcoxon; n=", nrow(paired), " pairs")
  } else {
    p <- suppressWarnings(wilcox.test(value ~ group, data = z, exact = FALSE)$p.value)
    method <- "patient-level Wilcoxon rank-sum"
  }
  c(p = p, method = method)
}

p_label <- function(p) {
  if (p < .001) "***  P<0.001" else if (p < .01) sprintf("**  P=%.3f", p) else if (p < .05) sprintf("*  P=%.3f", p) else sprintf("ns  P=%.3f", p)
}

message("[1/8] Fit the locked three-gene model and reconstruct cohorts")
expr <- read.csv(file.path(input_dir, "development_expression.csv"), row.names = 1, check.names = FALSE)
split <- read.csv(file.path(input_dir, "patient_split.csv"), check.names = FALSE)
dev <- as.data.frame(t(expr[genes, , drop = FALSE]), check.names = FALSE)
dev$sample_id <- rownames(dev)
dev <- merge(dev, split[, c("sample_id", "patient_id", "dataset", "group", "cohort")], by = "sample_id", sort = FALSE)
dev$group <- factor(ifelse(dev$group == "Treat", "PD", "NC"), levels = c("NC", "PD"))
dev$y <- as.integer(dev$group == "PD")
train <- dev[dev$cohort == "Training", , drop = FALSE]
test <- dev[dev$cohort == "Internal test", , drop = FALSE]
stopifnot(length(intersect(unique(train$patient_id), unique(test$patient_id))) == 0)
fit <- glm(y ~ AGT + CXCR4 + FOS, data = train, family = binomial())
train$probability <- predict(fit, newdata = train, type = "response")
test$probability <- predict(fit, newdata = test, type = "response")

ext_expr <- read.csv(file.path(input_dir, "GSE16134_expression.csv"), row.names = 1, check.names = FALSE)
ext_info <- read.csv(file.path(input_dir, "GSE16134_sample_info.csv"), check.names = FALSE)
ext <- as.data.frame(t(ext_expr[genes, , drop = FALSE]), check.names = FALSE)
ext$sample_id <- rownames(ext)
ext <- merge(ext, ext_info[, c("Sample", "Group", "Title")], by.x = "sample_id", by.y = "Sample", sort = FALSE)
ext$group <- factor(ifelse(ext$Group == "PD", "PD", "NC"), levels = c("NC", "PD"))
ext$y <- as.integer(ext$group == "PD")
ext$patient_id <- ifelse(grepl("patient [0-9]+", ext$Title, ignore.case = TRUE),
                         paste0("GSE16134_P", sub(".*patient ([0-9]+).*", "\\1", ext$Title, ignore.case = TRUE)), ext$sample_id)
# Use the Table S6 transport rule with the same training-fitted coefficients.
train_mu <- colMeans(train[, genes, drop = FALSE])
train_sd <- vapply(train[, genes, drop = FALSE], sd, numeric(1))
beta_raw <- coef(fit)[genes]
beta_z <- beta_raw * train_sd
intercept_z <- unname(coef(fit)[1] + sum(beta_raw * train_mu))
ext_z <- scale(as.matrix(ext[, genes, drop = FALSE]))
stopifnot(all(is.finite(ext_z)))
ext$probability <- plogis(intercept_z + drop(ext_z %*% beta_z))
stopifnot(abs(auc_rank(ext$y, ext$probability) - 0.9225449515905948) < 1e-10,
          abs(mean((ext$y - ext$probability)^2) - 0.14436794093654878) < 1e-10)
write.csv(ext[, c("sample_id", "patient_id", "group", "y", "probability")],
          file.path(out_dir, "GSE16134_fixed_predictions_S6.csv"), row.names = FALSE)
cohorts <- list(Training = train, `Internal test` = test, GSE16134 = ext)

metric_all <- read.csv(file.path(sensitivity_dir, "Table_AUC_cluster_bootstrap_CI.csv"), check.names = FALSE)
metric_three <- metric_all[metric_all$model == "Three_gene" & metric_all$cohort %in% c("Training", "Internal_test", "GSE16134"), ]
metric_three$Cohort <- sub("_", " ", metric_three$cohort)

message("[2/8] Figure 5A-B: exact additive SHAP values on the log-odds scale")
beta <- coef(fit)[genes]
center <- colMeans(train[, genes, drop = FALSE])
shap_wide <- sweep(as.matrix(train[, genes, drop = FALSE]), 2, center, "-")
shap_wide <- sweep(shap_wide, 2, beta, "*")
shap_long <- do.call(rbind, lapply(genes, function(g) {
  val <- train[[g]]; rng <- range(val)
  data.frame(Gene = g, SHAP = shap_wide[, g], Feature = (val - rng[1]) / diff(rng))
}))
importance <- aggregate(abs(SHAP) ~ Gene, shap_long, mean); names(importance)[2] <- "MeanAbs"
importance$Gene <- factor(importance$Gene, levels = importance$Gene[order(importance$MeanAbs)])
shap_long$Gene <- factor(shap_long$Gene, levels = levels(importance$Gene))

pA <- ggplot(shap_long, aes(SHAP, Gene, color = Feature)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = .5) +
  geom_point(position = position_jitter(width = 0, height = .12, seed = 20260902), size = 1.5, alpha = .82) +
  scale_color_gradientn(colors = c("#FDE725", "#E64B75", "#3B0F70"), limits = c(0, 1)) +
  labs(title = "Three-gene model: SHAP distribution", subtitle = "Exact additive contributions on the log-odds scale",
       x = "SHAP value (impact on log-odds)", y = NULL, color = "Feature value") + theme_pub(12)
pB <- ggplot(importance, aes(MeanAbs, Gene)) + geom_col(fill = "#12A7B5", width = .62) +
  labs(title = "Three-gene model: mean absolute SHAP", x = "Mean |SHAP value|", y = NULL) + theme_pub(12)
fA <- panel_png("Figure5_A_SHAP", pA, 5.7, 4.4)
fB <- panel_png("Figure5_B_SHAP_importance", pB, 5.0, 4.4)

message("[3/8] Figure 5C-E: ROC and expression panels")
roc_dev <- rbind(roc_df(train$y, train$probability, "Training"), roc_df(test$y, test$probability, "Internal test"))
roc_dev$Cohort <- factor(roc_dev$Cohort, levels = c("Training", "Internal test"))
lab_dev <- metric_three[metric_three$Cohort %in% c("Training", "Internal test"), ]
lab_dev$Cohort <- factor(lab_dev$Cohort, levels = levels(roc_dev$Cohort))
lab_dev$Label <- sprintf("AUC %.3f\n95%% CI %.3f-%.3f", lab_dev$AUC, lab_dev$CI_low, lab_dev$CI_high)
pC <- ggplot(roc_dev, aes(FPR, TPR, color = Cohort)) + geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_path(linewidth = 1.15) + geom_text(data = lab_dev, aes(.96, .08, label = Label), inherit.aes = FALSE, hjust = 1, size = 3.15) +
  facet_wrap(~Cohort, nrow = 1) + coord_equal() + scale_color_manual(values = cohort_cols, guide = "none") +
  labs(title = "Discrimination of the locked three-gene model", subtitle = "Patient-cluster bootstrap 95% confidence intervals",
       x = "1 - Specificity", y = "Sensitivity") + theme_pub(12)
fC <- panel_png("Figure5_C_internal_ROC", pC, 8.0, 4.4)

expr_long <- do.call(rbind, lapply(genes, function(g) data.frame(patient_id = train$patient_id, group = train$group, Gene = g, Expression = train[[g]])))
expr_long$Gene <- factor(expr_long$Gene, levels = genes)
expr_p <- do.call(rbind, lapply(genes, function(g) {
  d <- expr_long[expr_long$Gene == g, ]; z <- patient_test(d, "Expression")
  data.frame(Gene = g, P_value = as.numeric(z["p"]), Method = unname(z["method"]), Label = p_label(as.numeric(z["p"])), y = max(d$Expression) + .06 * diff(range(d$Expression)))
}))
pD <- ggplot(expr_long, aes(group, Expression, color = group)) +
  geom_jitter(width = .16, alpha = .67, size = 1.3) + geom_boxplot(width = .46, outlier.shape = NA, fill = NA, linewidth = .72) +
  geom_text(data = expr_p, aes(1.5, y, label = Label), inherit.aes = FALSE, size = 3.15) +
  facet_wrap(~Gene, nrow = 1, scales = "free_y") + scale_color_manual(values = group_cols) +
  labs(title = "Expression of the three retained biomarkers", subtitle = "Training cohort; patient-level paired tests where applicable",
       x = NULL, y = "Normalized expression", color = "Group") + theme_pub(12) + theme(legend.position = "right")
fD <- panel_png("Figure5_D_expression", pD, 8.8, 4.6)

pred223 <- read.csv(file.path(sensitivity_dir, "GSE223924_three_model_predictions.csv"), check.names = FALSE)
pred223$probability <- pred223$Three_gene_probability_MOR
roc_ext <- rbind(roc_df(ext$y, ext$probability, "GSE16134"), roc_df(pred223$y, pred223$probability, "GSE223924"))
roc_ext$Cohort <- factor(roc_ext$Cohort, levels = c("GSE16134", "GSE223924"))
m223 <- read.csv(file.path(sensitivity_dir, "Table_GSE223924_three_model_AUC.csv"), check.names = FALSE)
m223 <- m223[m223$model == "Three_gene" & m223$preprocessing == "Median_of_ratios_log2", ]
lab_ext <- rbind(
  data.frame(Cohort = "GSE16134", AUC = metric_three$AUC[metric_three$Cohort == "GSE16134"], CI_low = metric_three$CI_low[metric_three$Cohort == "GSE16134"], CI_high = metric_three$CI_high[metric_three$Cohort == "GSE16134"], n = nrow(ext), NC = sum(ext$y == 0), PD = sum(ext$y == 1)),
  data.frame(Cohort = "GSE223924", AUC = m223$AUC, CI_low = m223$CI_low, CI_high = m223$CI_high, n = nrow(pred223), NC = sum(pred223$y == 0), PD = sum(pred223$y == 1))
)
lab_ext$Cohort <- factor(lab_ext$Cohort, levels = levels(roc_ext$Cohort))
lab_ext$Label <- sprintf("n=%d; NC/PD=%d/%d\nAUC %.3f (95%% CI %.3f-%.3f)", lab_ext$n, lab_ext$NC, lab_ext$PD, lab_ext$AUC, lab_ext$CI_low, lab_ext$CI_high)
pE <- ggplot(roc_ext, aes(FPR, TPR, color = Cohort)) + geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_path(linewidth = 1.15) + geom_text(data = lab_ext, aes(.96, .07, label = Label), inherit.aes = FALSE, hjust = 1, size = 2.95) +
  facet_wrap(~Cohort, nrow = 1, labeller = as_labeller(c(GSE16134 = "GSE16134 (overlap; n=310)", GSE223924 = "GSE223924 (RNA-seq; n=20)"))) +
  coord_equal() + scale_color_manual(values = cohort_cols, guide = "none") +
  labs(title = "Additional cohort evaluation of the locked three-gene model", subtitle = "Full GSE16134 overlaps with development; RNA-seq evaluation is exploratory",
       x = "1 - Specificity", y = "Sensitivity") + theme_pub(12)
fE <- panel_png("Figure5_E_external_ROC", pE, 8.0, 4.6)

message("[4/8] Figure 5F-H: DCA, nomogram and calibration")
ths <- seq(.01, .80, .01)
single_pred <- list(`Three-gene model` = test$probability)
for (g in genes) {
  fg <- glm(as.formula(paste("y ~", g)), data = train, family = binomial())
  single_pred[[g]] <- predict(fg, newdata = test, type = "response")
}
dca_single <- do.call(rbind, lapply(names(single_pred), function(nm) {
  z <- net_benefit(test$y, single_pred[[nm]], ths)
  data.frame(Threshold = ths, NetBenefit = z$Model, Model = nm)
}))
ref_dca <- net_benefit(test$y, test$probability, ths)
pF <- ggplot(dca_single, aes(Threshold, NetBenefit, color = Model)) + geom_line(linewidth = .9) +
  geom_line(data = data.frame(Threshold = ths, NetBenefit = ref_dca$All), aes(Threshold, NetBenefit), inherit.aes = FALSE, color = "grey55", linetype = 2) +
  geom_hline(yintercept = 0, color = "black") + scale_color_manual(values = gene_cols) + coord_cartesian(ylim = c(-.12, 1.02)) +
  labs(title = "Three-gene model versus individual genes", subtitle = "Patient-separated internal test set",
       x = "Threshold probability", y = "Standardized net benefit", color = NULL) + theme_pub(11) + theme(legend.position = c(.78, .75), legend.background = element_rect(fill = "white", color = "grey70"))
fF <- panel_png("Figure5_F_single_vs_combined_DCA", pF, 5.9, 4.8)

## Graphical nomogram: coefficient-weighted gene points sum to a calibrated total-points/probability axis.
gene_quantiles <- lapply(genes, function(g) as.numeric(quantile(train[[g]], seq(0, 1, length.out = 6))))
names(gene_quantiles) <- genes
gene_contrib <- lapply(genes, function(g) beta[g] * gene_quantiles[[g]]); names(gene_contrib) <- genes
gene_min <- vapply(gene_contrib, min, numeric(1))
gene_range <- vapply(gene_contrib, function(x) diff(range(x)), numeric(1))
point_scale <- 100 / max(gene_range)
nom_rows <- do.call(rbind, lapply(seq_along(genes), function(i) {
  g <- genes[i]
  data.frame(Gene=g, Expression=gene_quantiles[[g]], Points=(gene_contrib[[g]]-gene_min[g])*point_scale, y=length(genes)-i+2)
}))
max_total <- sum(gene_range * point_scale)
base_lp <- unname(coef(fit)[1] + sum(gene_min))
total_ticks <- seq(0, max_total, length.out=6)
prob_vals <- c(.05,.10,.25,.50,.75,.90,.95)
prob_points <- (qlogis(prob_vals)-base_lp)*point_scale
keep_prob <- prob_points >= 0 & prob_points <= max_total
prob_axis <- data.frame(Points=prob_points[keep_prob],Probability=prob_vals[keep_prob],y=.65)
nom_rows$Gene <- factor(nom_rows$Gene, levels = rev(genes))
pG <- ggplot(nom_rows, aes(Points, y)) +
  geom_segment(data = do.call(rbind, lapply(split(nom_rows, nom_rows$Gene), function(z) data.frame(x=min(z$Points), xend=max(z$Points), y=z$y[1], yend=z$y[1]))), aes(x=x,xend=xend,y=y,yend=yend), inherit.aes=FALSE, linewidth=.65) +
  geom_point(size = 1.7, color = "#4056C9") +
  geom_text(aes(label = sprintf("%.1f", Expression)), vjust = -1.0, size = 2.5) +
  annotate("segment",x=0,xend=max_total,y=1.55,yend=1.55,linewidth=.75) +
  annotate("point",x=total_ticks,y=1.55,size=1.5) +
  annotate("text",x=total_ticks,y=1.55,label=sprintf("%.0f",total_ticks),vjust=-1,size=2.45) +
  annotate("text",x=-.025*max_total,y=1.55,label="Total points",hjust=1,fontface="bold",size=3.0) +
  geom_segment(data=prob_axis,aes(x=Points,xend=Points,y=.58,yend=.72),inherit.aes=FALSE,linewidth=.5) +
  geom_text(data=prob_axis,aes(Points,.45,label=sprintf("%.2f",Probability)),inherit.aes=FALSE,size=2.45) +
  annotate("segment",x=0,xend=max_total,y=.65,yend=.65,linewidth=.75) +
  annotate("text",x=-.025*max_total,y=.65,label="PD probability",hjust=1,fontface="bold",size=3.0) +
  scale_y_continuous(breaks=seq(2,length(genes)+1),labels=rev(genes),limits=c(.25,length(genes)+1.55)) +
  scale_x_continuous(limits = c(-.12*max_total, 1.02*max_total), breaks = seq(0,max_total,length.out=6), labels=function(x)sprintf("%.0f",x), position = "top") +
  labs(title = "Three-gene nomogram", subtitle = "Expression values mapped to points, total points and PD probability", x = "Points", y = NULL) + theme_pub(10.5) +
  theme(panel.grid.major.x = element_line(color = "grey90", linewidth = .35), axis.line.y = element_blank(), axis.ticks.y = element_blank())
fG <- panel_png("Figure5_G_nomogram", pG, 6.6, 4.8)

cal_dat <- rbind(calibration_bins(train, "Training"), calibration_bins(test, "Internal test"))
cal_dat$Cohort <- factor(cal_dat$Cohort, levels = c("Training", "Internal test"))
cal_lab <- do.call(rbind, lapply(c("Training", "Internal test"), function(nm) {
  d <- if (nm == "Training") train else test; cs <- calibration_stats(d$y, d$probability)
  data.frame(Cohort = nm, Label = sprintf("Brier %.3f\nIntercept %.2f; slope %.2f", cs[1], cs[2], cs[3]))
}))
cal_lab$Cohort <- factor(cal_lab$Cohort, levels = levels(cal_dat$Cohort))
pH <- ggplot(cal_dat, aes(probability, y, color = Cohort)) + geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
  geom_line(linewidth = 1) + geom_point(size = 2.2) + geom_text(data = cal_lab, aes(.04, .96, label = Label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.75) +
  facet_wrap(~Cohort, nrow = 1) + coord_equal(xlim = c(0,1), ylim = c(0,1)) + scale_color_manual(values = cohort_cols, guide = "none") +
  labs(title = "Calibration of locked predictions", subtitle = "Observed event rate versus predicted probability in quantile bins",
       x = "Mean predicted probability", y = "Observed periodontitis proportion") + theme_pub(11)
fH <- panel_png("Figure5_H_calibration", pH, 7.2, 4.8)

message("[5/8] Figure S2D-G: selection audit and locked-model sensitivity panels")
screen_auc <- data.frame(
  Model=c("LASSO","SVM","NNET","RF","KNN","GBM","DT","GLM"),
  AUC=c(.8702791461,.8653530378,.8538587849,.8357963875,.8218390805,.8193760263,.7750410509,.7110016420)
)
screen_auc$Model <- factor(screen_auc$Model, levels=rev(screen_auc$Model[order(screen_auc$AUC,decreasing=TRUE)]))
pC2 <- ggplot(screen_auc,aes(AUC,Model)) +
  geom_segment(aes(x=.5,xend=AUC,y=Model,yend=Model),color="grey75",linewidth=1.0) +
  geom_point(aes(color=AUC),size=4.2) +
  geom_text(aes(label=sprintf("%.3f",AUC)),hjust=-.35,size=3.4) +
  geom_vline(xintercept=.5,linetype=2,color="grey55") +
  scale_color_gradient(low="#86C7C4",high="#235AA6",guide="none") +
  scale_x_continuous(limits=c(.48,.91),breaks=seq(.5,.9,.1)) +
  labs(title="Eight-algorithm discrimination",subtitle="Locked independent test set (n = 50)",x="Area under the ROC curve",y=NULL) + theme_pub(12)
fC2 <- panel_png("FigureS2_C_test_AUC_ranking",pC2,5.6,4.7)

pD2 <- ggplot() +
  annotate("rect", xmin=.04,xmax=.96,ymin=.70,ymax=.96,fill="#EEF4FB",color="#5A7DA5",linewidth=.7) +
  annotate("text", x=.5,y=.90,label="Four-method top-20 consensus",fontface="bold",size=4.2) +
  annotate("text", x=.5,y=.80,label="AGT   CXCR4   FOS   EPHX2",size=4.2) +
  annotate("segment", x=.5,xend=.5,y=.69,yend=.57,arrow=arrow(length=grid::unit(.14,"in")),linewidth=.7) +
  annotate("rect", xmin=.04,xmax=.96,ymin=.30,ymax=.57,fill="#FFF7E8",color="#CC8B2C",linewidth=.7) +
  annotate("text", x=.5,y=.51,label="Training-only parsimony check",fontface="bold",size=4.2) +
  annotate("text", x=.5,y=.41,label="Incremental EPHX2: LRT P = 0.595\nAIC: 3-gene 94.999 vs 4-gene 96.717",size=3.6) +
  annotate("segment", x=.5,xend=.5,y=.29,yend=.18,arrow=arrow(length=grid::unit(.14,"in")),linewidth=.7) +
  annotate("rect", xmin=.04,xmax=.96,ymin=.02,ymax=.18,fill="#EAF7EF",color="#2C8B57",linewidth=.8) +
  annotate("text", x=.5,y=.10,label="Final locked model: AGT + CXCR4 + FOS",fontface="bold",size=4.1) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="off") + labs(title="Feature-selection audit trail") + theme_void(base_size=12) + theme(plot.title=element_text(face="bold",hjust=.5))
fD2 <- panel_png("FigureS2_D_selection_audit", pD2, 6.0, 5.0)

prob_long <- do.call(rbind, lapply(names(cohorts), function(nm) {
  d <- cohorts[[nm]]; data.frame(Cohort=nm, patient_id=d$patient_id, group=d$group, probability=d$probability)
}))
prob_long$Cohort <- factor(prob_long$Cohort, levels = names(cohorts))
prob_p <- do.call(rbind, lapply(names(cohorts), function(nm) {
  d <- prob_long[prob_long$Cohort == nm, ]; z <- patient_test(d, "probability")
  data.frame(Cohort=nm, Label=p_label(as.numeric(z["p"])), P_value=as.numeric(z["p"]), Method=unname(z["method"]), y=1.025)
}))
prob_p$Cohort <- factor(prob_p$Cohort, levels=names(cohorts))
pE2 <- ggplot(prob_long, aes(group, probability, color=group)) + geom_jitter(width=.15,alpha=.65,size=1.2) +
  geom_boxplot(width=.45,outlier.shape=NA,fill=NA,linewidth=.7) + geom_text(data=prob_p,aes(1.5,y,label=Label),inherit.aes=FALSE,size=3.0) +
  facet_wrap(~Cohort,nrow=1) + scale_color_manual(values=group_cols) + scale_y_continuous(limits=c(0,1.08),breaks=seq(0,1,.2)) +
  labs(title="Locked three-gene model probabilities",subtitle="Patient-level paired Wilcoxon tests when paired sites are available",x=NULL,y="Predicted probability",color="Group") + theme_pub(11)
fE2 <- panel_png("FigureS2_E_probabilities", pE2, 9.6, 5.0)

heat <- data.frame(
  Gene=rep(c("CXCR4","FOS"),each=7),
  Module=rep(c("MEblue","MEbrown","MEgreen","MEgrey","MEred","MEturquoise","MEyellow"),2),
  r=c(.87,.12,.22,.01,.60,-.71,-.15, .39,-.10,.22,.16,.41,-.42,.18)
)
heat$Gene <- factor(heat$Gene,levels=rev(c("CXCR4","FOS")))
heat$Module <- factor(heat$Module,levels=c("MEblue","MEbrown","MEgreen","MEgrey","MEred","MEturquoise","MEyellow"))
pF2 <- ggplot(heat,aes(Module,Gene,fill=r)) + geom_tile(color="grey70",linewidth=.45) + geom_text(aes(label=sprintf("%.2f",r)),size=3.3) +
  scale_fill_gradient2(low="#3973B7",mid="#FFF9C4",high="#EF3B2C",midpoint=0,limits=c(-.9,.9),name="r") +
  labs(title="Biomarker-module eigengene correlations",subtitle="Original WGCNA correlations; C3 row removed with model revision",x=NULL,y=NULL) +
  theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5,size=9,color="grey35"),axis.text.x=element_text(angle=90,vjust=.5,hjust=1,color="black"),axis.text.y=element_text(color="black"),panel.grid=element_blank())
fF2 <- panel_png("FigureS2_F_WGCNA_heatmap", pF2, 6.0, 4.8)

dca_all <- do.call(rbind,lapply(names(cohorts),function(nm){z<-net_benefit(cohorts[[nm]]$y,cohorts[[nm]]$probability,ths);z$Cohort<-nm;z}))
dca_long <- reshape(dca_all,varying=c("Model","All","None"),v.names="NetBenefit",timevar="Strategy",times=c("Three-gene model","Treat all","Treat none"),direction="long")
dca_long$Cohort <- factor(dca_long$Cohort,levels=names(cohorts))
dca_long$Strategy <- factor(dca_long$Strategy,levels=c("Three-gene model","Treat all","Treat none"))
pG2 <- ggplot(dca_long,aes(Threshold,NetBenefit,color=Strategy,linetype=Strategy)) + geom_line(linewidth=.9) + facet_wrap(~Cohort,nrow=1,scales="free_y") +
  scale_color_manual(values=c(`Three-gene model`=model_col,`Treat all`="grey50",`Treat none`="black")) +
  scale_linetype_manual(values=c(`Three-gene model`=1,`Treat all`=2,`Treat none`=1)) +
  labs(title="Decision curve analysis",subtitle="Locked predictions; standardized net benefit",x="Threshold probability",y="Standardized net benefit",color=NULL,linetype=NULL) +
  theme_pub(11) + theme(legend.position="bottom")
fG2 <- panel_png("FigureS2_G_DCA", pG2, 9.6, 4.8)

message("[6/8] Compose publication figures")
add_label <- function(img, lab, size=82) {
  image_annotate(img, lab, gravity="northwest", location="+18+8", size=size, font="Arial", weight=700, color="black")
}
fit_box <- function(path, width, height, lab=NULL) {
  x <- if (inherits(path, "magick-image")) path else image_read(path)
  x <- image_trim(x, fuzz=2)
  x <- image_resize(x, paste0(width,"x",height,">"))
  x <- image_extent(x, paste0(width,"x",height), gravity="center", color="white")
  if (!is.null(lab)) x <- add_label(x, lab, max(52, round(height*.08)))
  x
}
row_join <- function(items) image_append(do.call(c, items), stack=FALSE)
stack_join <- function(items) image_append(do.call(c, items), stack=TRUE)

fig5 <- stack_join(list(
  row_join(list(fit_box(fA,1900,1600,"A"),fit_box(fB,1500,1600,"B"),fit_box(fC,2600,1600,"C"))),
  row_join(list(fit_box(fD,3300,1700,"D"),fit_box(fE,2700,1700,"E"))),
  row_join(list(fit_box(fF,1900,1700,"F"),fit_box(fG,2100,1700,"G"),fit_box(fH,2000,1700,"H")))
))
fig5 <- image_extent(fig5,"6000x5000",gravity="center",color="white")

oldA <- image_read("source_S2_A.png")
oldB <- image_read("source_S2_B.png")
figs2 <- stack_join(list(
  row_join(list(fit_box(oldA,1900,1550,"A"),fit_box(oldB,1900,1550,"B"),fit_box(fC2,2200,1550,"C"))),
  row_join(list(fit_box(fD2,2200,1650,"D"),fit_box(fE2,3800,1650,"E"))),
  row_join(list(fit_box(fF2,2200,1650,"F"),fit_box(fG2,3800,1650,"G")))
))
figs2 <- image_extent(figs2,"6000x4850",gravity="center",color="white")

write_outputs <- function(img, stem) {
  image_write(img, file.path(out_dir,paste0(stem,".png")), format="png")
  image_write(img, file.path(out_dir,paste0(stem,".tif")), format="tiff", compression="lzw")
  image_write(img, file.path(out_dir,paste0(stem,".pdf")), format="pdf")
  image_write(img, file.path(figure_dir,paste0(stem,".png")), format="png")
  image_write(img, file.path(figure_dir,paste0(stem,".tif")), format="tiff", compression="lzw")
}
write_outputs(fig5,"Figure5_three_gene_20260902")
write_outputs(figs2,"FigureS2_three_gene_20260902")

message("[7/8] Save numerical outputs and reproducibility records")
write.csv(data.frame(Term=names(coef(fit)),Coefficient=unname(coef(fit)),Odds_ratio=exp(unname(coef(fit)))),file.path(out_dir,"Table_three_gene_coefficients.csv"),row.names=FALSE)
write.csv(metric_three,file.path(out_dir,"Table_three_gene_AUC_cluster_bootstrap_CI.csv"),row.names=FALSE)
write.csv(expr_p,file.path(out_dir,"Table_expression_patient_level_tests.csv"),row.names=FALSE)
write.csv(prob_p,file.path(out_dir,"Table_probability_patient_level_tests.csv"),row.names=FALSE)
write.csv(dca_long,file.path(out_dir,"Table_three_gene_DCA_values.csv"),row.names=FALSE)
write.csv(heat,file.path(out_dir,"Table_biomarker_module_correlations.csv"),row.names=FALSE)
saveRDS(fit,file.path(out_dir,"Locked_three_gene_logistic_model.rds"))
writeLines(c(
  "FIGURE 5 AND FIGURE S2 THREE-GENE REVISION",
  "Final locked genes: AGT, CXCR4, FOS.",
  "EPHX2 was not retained because its incremental training-set likelihood-ratio P value was 0.595 and the three-gene model had lower AIC (94.999 versus 96.717).",
  "Figure S2A-B retain the original eight-algorithm residual diagnostics; Figure S2C reports the formal fixed-seed 41-gene rerun's independent-test AUC ranking. Downstream locked-model panels D-G were regenerated.",
  "SHAP values in Figure 5A-B are exact centered additive contributions on the logistic link (log-odds) scale.",
  "No patient appears in both training and internal test sets. External predictions were generated without model refitting or threshold re-optimization.",
  "GSE223924 uses the prespecified median-of-ratios log2 preprocessing and is presented as exploratory cross-platform validation."
),file.path(out_dir,"README_revision_scope.txt"))
capture.output(sessionInfo(),file=file.path(out_dir,"R_sessionInfo.txt"))
pk <- data.frame(Package=required_packages,Version=vapply(required_packages,function(p)as.character(packageVersion(p)),character(1)))
write.csv(pk,file.path(out_dir,"R_package_versions.csv"),row.names=FALSE)

message("[8/8] Completed: ", out_dir)
