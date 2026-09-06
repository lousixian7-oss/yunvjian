if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent=TRUE)
options(stringsAsFactors=FALSE)
set.seed(20260905)
user_library <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
.libPaths(c(user_library, "G:/gurobi/Rlib45", .libPaths()))
stopifnot(requireNamespace("ggplot2", quietly=TRUE))
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2L) args <- c("G:/1Yunvjian/1A\u7389\u5973\u714e/code", "G:/1Yunvjian/\u673a\u5668\u5b66\u4e60/8.10/\u7ed3\u679c/tripod_20260905/probability_reconciliation")
archive <- normalizePath(args[1], winslash="/", mustWork=TRUE)
out_dir <- normalizePath(args[2], winslash="/", mustWork=TRUE)
dir.create(file.path(out_dir,"panels"),showWarnings=FALSE)
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


expr <- read.csv(gzfile(file.path(archive,"coredata/bulk/development_expression.csv.gz")),row.names=1,check.names=FALSE)
split <- read.csv(file.path(archive,"coredata/bulk/patient_split.csv"))
dev <- as.data.frame(t(expr[genes,,drop=FALSE])); dev$sample_id <- rownames(dev)
dev <- merge(dev,split,by="sample_id",sort=FALSE)
dev$group <- factor(ifelse(dev$group=="Treat","PD","NC"),levels=c("NC","PD"));dev$y <- as.integer(dev$group=="PD")
train <- dev[dev$cohort=="Training",];test <- dev[dev$cohort=="Internal test",]
stopifnot(nrow(train)==120L,nrow(test)==52L,length(intersect(train$patient_id,test$patient_id))==0L)
fit <- readRDS(file.path(archive,"coredata/model/Locked_three_gene_logistic_model.rds"))
stopifnot(max(abs(coef(fit)-c(-28.1615557900581,2.02555945620692,1.70618768640773,0.412233157693864)))<1e-10)
train$probability <- predict(fit,newdata=train,type="response")
test$probability <- predict(fit,newdata=test,type="response")
ex <- read.csv(gzfile(file.path(archive,"coredata/bulk/GSE16134_expression.csv.gz")),row.names=1,check.names=FALSE)
ext <- as.data.frame(t(ex[genes,,drop=FALSE]));ext$sample_id<-rownames(ext)
info<-read.csv(file.path(archive,"coredata/bulk/GSE16134_sample_info.csv"))
ext<-merge(ext,info[,c("Sample","Group","Title")],by.x="sample_id",by.y="Sample",sort=FALSE)
ext$group<-factor(ext$Group,levels=c("NC","PD"));ext$y<-as.integer(ext$group=="PD")
ext$patient_id<-paste0("GSE16134_P",sub(".*patient ([0-9]+).*","\\1",ext$Title,ignore.case=TRUE))
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
cohorts<-list(Training=train,`Internal test`=test,GSE16134=ext)
ci<-read.csv(file.path(archive,"coredata/model/corrected_release/Table_three_gene_AUC_cluster_bootstrap_CI.csv"))
if(!"Cohort" %in% names(ci)) ci$Cohort<-sub("_"," ",ci$cohort)
metric_three<-ci
sensitivity_dir<-file.path(archive,"coredata/model/fixed_model")
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

ths<-seq(.01,.80,.01)
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


all_pred <- do.call(rbind,lapply(names(cohorts),function(nm){d<-cohorts[[nm]];data.frame(Cohort=nm,sample_id=d$sample_id,patient_id=d$patient_id,group=d$group,y=d$y,probability=d$probability)}))
all_pred<-rbind(all_pred,data.frame(Cohort="GSE223924",sample_id=pred223$sample_id,patient_id=pred223$sample_id,group=pred223$group,y=pred223$y,probability=pred223$probability))
write.csv(all_pred,file.path(out_dir,"Unified_fixed_predictions.csv"),row.names=FALSE)
stats<-do.call(rbind,lapply(split(all_pred,all_pred$Cohort),function(d){v<-calibration_stats(d$y,d$probability);data.frame(Cohort=d$Cohort[1],AUC=auc_rank(d$y,d$probability),Brier=unname(v[1]),Calibration_intercept=unname(v[2]),Calibration_slope=unname(v[3]))}))
expected<-read.csv(file.path(out_dir,"S6_expected_metrics.csv"),check.names=FALSE)
for(i in seq_len(nrow(expected))){r<-stats[stats$Cohort==expected$Cohort[i],];stopifnot(nrow(r)==1L,max(abs(as.numeric(r[1,-1])-as.numeric(expected[i,-1])))<1e-8)}
write.csv(stats,file.path(out_dir,"Verified_S6_metrics.csv"),row.names=FALSE)
write.csv(prob_p,file.path(out_dir,"Table_probability_patient_level_tests.csv"),row.names=FALSE)
write.csv(dca_long,file.path(out_dir,"Table_three_gene_DCA_values.csv"),row.names=FALSE)
write.csv(data.frame(Gene=genes,Training_mean=train_mu,Training_SD=train_sd,External_mean=attr(ext_z,"scaled:center"),External_SD=attr(ext_z,"scaled:scale"),Standardized_coefficient=beta_z),file.path(out_dir,"External_scaling_constants.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out_dir,"R_sessionInfo.txt"))
write.csv(data.frame(package="ggplot2",version=as.character(packageVersion("ggplot2"))),file.path(out_dir,"R_package_versions.csv"),row.names=FALSE)
print(stats);message("All AUC, Brier and calibration results match Table S6. No model refitting or recalibration performed.")
