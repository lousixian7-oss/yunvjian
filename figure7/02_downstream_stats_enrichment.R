user_lib <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

set.seed(12345)
output_dir <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final/hdWGCNA_rerun_20260824"
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
object_file <- file.path(output_dir, "hdWGCNA_object_final.rds")
old_hub_file <- file.path(output_dir, "reference_inputs", "old_Fibroblasts1_hubgenes.csv")
gmt_legacy_file <- file.path(output_dir, "reference_inputs", "kegg_legacy_2026_1.gmt")
gmt_medicus_file <- file.path(output_dir, "reference_inputs", "kegg_medicus_2026_1.gmt")

log_con <- file(file.path(output_dir, "02_downstream_stats_enrichment.log"), open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Start:", format(Sys.time()), "\n")
obj <- readRDS(object_file)
DefaultAssay(obj) <- "RNA"

mes <- GetMEs(obj, harmonized = TRUE)
mods <- setdiff(colnames(mes), "grey")
stopifnot(identical(rownames(mes), colnames(obj)))

expr <- FetchData(obj, vars = c("AGT", "PI16"), layer = "data")
meta <- obj[[]][, c("orig.ident", "sample", "dataset", "group", "state_final"), drop = FALSE]
cell_df <- cbind(
  cell = rownames(meta),
  meta,
  AGT = expr[rownames(meta), "AGT"],
  PI16 = expr[rownames(meta), "PI16"],
  mes[rownames(meta), mods, drop = FALSE]
) %>% as.data.frame(check.names = FALSE)

module_summary <- cell_df %>%
  group_by(sample, dataset, group) %>%
  summarise(
    n_fibro = n(),
    AGT_mean = mean(AGT),
    AGT_positive_fraction = mean(AGT > 0),
    PI16_mean = mean(PI16),
    PI16_positive_fraction = mean(PI16 > 0),
    across(all_of(mods), mean),
    .groups = "drop"
  )

write.csv(module_summary, file.path(table_dir, "05_sample_level_AGT_PI16_module_hMEs.csv"), row.names = FALSE)
eligible <- module_summary %>% filter(n_fibro >= 50)

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

cor_rows <- lapply(mods, function(mod) {
  a <- safe_cor(eligible[[mod]], eligible$AGT_mean)
  b <- safe_cor(eligible[[mod]], eligible$AGT_positive_fraction)
  data.frame(
    module = mod,
    rho_AGT_mean = a[["rho"]], p_AGT_mean = a[["p"]], n_AGT_mean = a[["n"]],
    rho_AGT_positive_fraction = b[["rho"]], p_AGT_positive_fraction = b[["p"]],
    n_AGT_positive_fraction = b[["n"]]
  )
})
cor_table <- bind_rows(cor_rows) %>%
  mutate(
    p_adjust_AGT_mean = p.adjust(p_AGT_mean, method = "BH"),
    p_adjust_AGT_positive_fraction = p.adjust(p_AGT_positive_fraction, method = "BH")
  ) %>%
  arrange(desc(rho_AGT_mean))
write.csv(cor_table, file.path(table_dir, "06_AGT_module_spearman_all_modules.csv"), row.names = FALSE)

selected_module <- cor_table$module[which.max(cor_table$rho_AGT_mean)]
cat("Selected AGT-associated module:", selected_module, "\n")
print(cor_table)

state_summary <- cell_df %>%
  group_by(sample, dataset, group, state_final) %>%
  summarise(n_cells = n(), across(all_of(mods), mean), .groups = "drop") %>%
  filter(n_cells >= 10)
write.csv(state_summary, file.path(table_dir, "07_sample_state_module_hMEs_n10.csv"), row.names = FALSE)

state_global <- bind_rows(lapply(mods, function(mod) {
  d <- state_summary %>% filter(is.finite(.data[[mod]]))
  kw <- kruskal.test(d[[mod]] ~ d$state_final)
  data.frame(module = mod, statistic = unname(kw$statistic), p = kw$p.value, n_sample_state = nrow(d))
})) %>% mutate(p_adjust = p.adjust(p, method = "BH"))
write.csv(state_global, file.path(table_dir, "08_module_state_kruskal_tests.csv"), row.names = FALSE)

state_pairs <- combn(sort(unique(as.character(state_summary$state_final))), 2, simplify = FALSE)
pairwise_state <- bind_rows(lapply(mods, function(mod) {
  bind_rows(lapply(state_pairs, function(pair) {
    d <- state_summary %>% filter(state_final %in% pair)
    wt <- suppressWarnings(wilcox.test(d[[mod]] ~ d$state_final, exact = FALSE))
    data.frame(module = mod, state1 = pair[1], state2 = pair[2], W = unname(wt$statistic), p = wt$p.value)
  }))
})) %>% group_by(module) %>% mutate(p_adjust = p.adjust(p, method = "BH")) %>% ungroup()
write.csv(pairwise_state, file.path(table_dir, "09_module_state_pairwise_wilcoxon_BH.csv"), row.names = FALSE)

group_test <- suppressWarnings(wilcox.test(eligible[[selected_module]] ~ eligible$group, exact = FALSE))
group_lm_data <- eligible %>%
  mutate(score_z = as.numeric(scale(.data[[selected_module]])), group = factor(group, levels = c("NC", "PD")))
group_lm <- lm(score_z ~ dataset + group, data = group_lm_data)
group_lm_coef <- as.data.frame(summary(group_lm)$coefficients)
group_lm_coef$term <- rownames(group_lm_coef)
rownames(group_lm_coef) <- NULL
write.csv(group_lm_coef, file.path(table_dir, "10_selected_module_dataset_adjusted_group_model.csv"), row.names = FALSE)
write.csv(data.frame(
  module = selected_module,
  wilcoxon_W = unname(group_test$statistic),
  wilcoxon_p = group_test$p.value,
  n_NC = sum(eligible$group == "NC"), n_PD = sum(eligible$group == "PD")
), file.path(table_dir, "11_selected_module_PD_NC_wilcoxon.csv"), row.names = FALSE)

hub <- read.csv(file.path(table_dir, "03_hub_genes_top100_each_module.csv"), check.names = FALSE)
new_selected_hubs <- hub %>% filter(module == selected_module) %>% arrange(desc(kME))
write.csv(new_selected_hubs, file.path(table_dir, "12_selected_AGT_module_hub_genes.csv"), row.names = FALSE)

old_hubs <- read.csv(old_hub_file, check.names = FALSE)$gene_name
overlap <- hub %>% group_by(module) %>% summarise(
  new_top_n = n(), old_top_n = length(old_hubs), overlap_n = sum(gene_name %in% old_hubs),
  jaccard = overlap_n / (new_top_n + old_top_n - overlap_n),
  overlap_genes = paste(sort(intersect(gene_name, old_hubs)), collapse = ";"), .groups = "drop"
) %>% arrange(desc(jaccard))
write.csv(overlap, file.path(table_dir, "13_old_Fibroblasts1_vs_new_modules_top_hub_overlap.csv"), row.names = FALSE)

technical_pattern <- "^(MT-|RPS|RPL|HSP|MALAT1$|FOS$|FOSB$|JUN$|JUNB$|JUND$|MKI67$|TOP2A$)"
independent_genes <- new_selected_hubs %>%
  filter(!gene_name %in% c("AGT", "PI16"), !grepl(technical_pattern, gene_name, ignore.case = TRUE)) %>%
  slice_head(n = 25) %>% pull(gene_name)
if (length(independent_genes) < 10) {
  independent_genes <- new_selected_hubs %>%
    filter(!gene_name %in% c("AGT", "PI16")) %>% slice_head(n = 25) %>% pull(gene_name)
}
write.csv(data.frame(module = selected_module, rank = seq_along(independent_genes), gene = independent_genes),
          file.path(table_dir, "14_independent_module_score_genes_excluding_AGT_PI16_technical.csv"), row.names = FALSE)

obj <- AddModuleScore(obj, features = list(independent_genes), name = "AGT_independent_module", seed = 12345)
safe_score_col <- "AGT_independent_module1"
score_df <- obj[[]] %>%
  mutate(cell = rownames(.), state_final = as.character(state_final)) %>%
  dplyr::select(cell, sample, dataset, group, state_final, all_of(safe_score_col)) %>%
  group_by(sample, dataset, group, state_final) %>%
  summarise(n_cells = n(), independent_module_score = mean(.data[[safe_score_col]]), .groups = "drop") %>%
  filter(n_cells >= 10)
write.csv(score_df, file.path(table_dir, "15_independent_module_score_by_sample_state_n10.csv"), row.names = FALSE)

score_kw <- kruskal.test(independent_module_score ~ state_final, data = score_df)
score_pairs <- bind_rows(lapply(state_pairs, function(pair) {
  d <- score_df %>% filter(state_final %in% pair)
  wt <- suppressWarnings(wilcox.test(independent_module_score ~ state_final, data = d, exact = FALSE))
  data.frame(state1 = pair[1], state2 = pair[2], W = unname(wt$statistic), p = wt$p.value)
})) %>% mutate(p_adjust = p.adjust(p, method = "BH"))
write.csv(score_pairs, file.path(table_dir, "16_independent_score_state_pairwise_wilcoxon_BH.csv"), row.names = FALSE)
write.csv(data.frame(statistic = unname(score_kw$statistic), p = score_kw$p.value),
          file.path(table_dir, "17_independent_score_state_kruskal.csv"), row.names = FALSE)

state_levels <- c("PI16_Fib", "Adventitial_Fib", "ECM_Fib", "Activated_Fib", "Inflammatory_Fib")
module_long <- state_summary %>%
  dplyr::select(sample, state_final, all_of(mods)) %>%
  pivot_longer(all_of(mods), names_to = "module", values_to = "hME") %>%
  group_by(module) %>% mutate(hME_z = as.numeric(scale(hME))) %>% ungroup() %>%
  group_by(module, state_final) %>% summarise(mean_hME_z = mean(hME_z), .groups = "drop") %>%
  left_join(state_global %>% dplyr::select(module, p_adjust), by = "module") %>%
  mutate(
    state_final = factor(state_final, levels = state_levels),
    module = factor(module, levels = rev(mods)),
    sig = case_when(p_adjust < 0.001 ~ "***", p_adjust < 0.01 ~ "**", p_adjust < 0.05 ~ "*", TRUE ~ "ns")
  )

p_d <- ggplot(module_long, aes(state_final, module, fill = mean_hME_z)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = module_long %>% group_by(module) %>% slice(1),
            aes(x = Inf, y = module, label = sig), inherit.aes = FALSE, hjust = 1.2, size = 4) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "mean z-scored hME",
       title = "Module activity across revised fibroblast states",
       subtitle = "GSM×state means (≥10 cells); stars denote BH-adjusted Kruskal–Wallis tests") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
ggsave(file.path(figure_dir, "Figure7D_sample_level_module_state_heatmap.pdf"), p_d, width = 9, height = 5.5, bg = "white")
ggsave(file.path(figure_dir, "Figure7D_sample_level_module_state_heatmap.png"), p_d, width = 9, height = 5.5, dpi = 400, bg = "white")

sel_cor <- cor_table %>% filter(module == selected_module)
p_e <- ggplot(eligible, aes(x = .data[[selected_module]], y = AGT_mean, color = dataset, shape = group)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "grey35", linewidth = 0.7) +
  labs(x = paste0(selected_module, " mean harmonized eigengene"), y = "Mean AGT expression",
       title = "Sample-level association between AGT and the selected module",
       subtitle = sprintf("Spearman rho = %.2f, p = %.3g, BH-adjusted p = %.3g; GSMs with ≥50 fibroblasts",
                          sel_cor$rho_AGT_mean, sel_cor$p_AGT_mean, sel_cor$p_adjust_AGT_mean)) +
  theme_bw(base_size = 11)
ggsave(file.path(figure_dir, "Figure7E_AGT_selected_module_sample_correlation.pdf"), p_e, width = 7, height = 5, bg = "white")
ggsave(file.path(figure_dir, "Figure7E_AGT_selected_module_sample_correlation.png"), p_e, width = 7, height = 5, dpi = 400, bg = "white")

p_f <- ggplot(score_df, aes(x = factor(state_final, levels = state_levels), y = independent_module_score, fill = state_final)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = dataset), width = 0.14, size = 2, alpha = 0.85) +
  labs(x = NULL, y = "Independent module score", fill = "state", color = "dataset",
       title = paste0(selected_module, " program score across fibroblast states"),
       subtitle = sprintf("Top 25 hubs excluding AGT, PI16 and technical genes; Kruskal–Wallis p = %.3g", score_kw$p.value)) +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(figure_dir, "Figure7F_independent_module_score_by_state.pdf"), p_f, width = 8.5, height = 5.5, bg = "white")
ggsave(file.path(figure_dir, "Figure7F_independent_module_score_by_state.png"), p_f, width = 8.5, height = 5.5, dpi = 400, bg = "white")

lm_group_p <- group_lm_coef %>% filter(term == "groupPD") %>% pull(`Pr(>|t|)`)
p_g <- ggplot(eligible, aes(group, .data[[selected_module]], fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(aes(color = dataset), width = 0.12, size = 2.4, alpha = 0.9) +
  scale_fill_manual(values = c(NC = "#36BFC4", PD = "#F8766D")) +
  labs(x = NULL, y = paste0(selected_module, " mean harmonized eigengene"),
       title = paste0(selected_module, " activity in NC and PD"),
       subtitle = sprintf("Sample-level Wilcoxon p = %.3g; dataset-adjusted linear-model p = %.3g; ≥50 fibroblasts/GSM",
                          group_test$p.value, lm_group_p), color = "dataset", fill = "group") +
  theme_bw(base_size = 11)
ggsave(file.path(figure_dir, "Figure7G_selected_module_PD_NC_sample_level.pdf"), p_g, width = 6.5, height = 5, bg = "white")
ggsave(file.path(figure_dir, "Figure7G_selected_module_PD_NC_sample_level.png"), p_g, width = 6.5, height = 5, dpi = 400, bg = "white")

# Enrichment of the top 100 hub genes from the newly selected AGT-associated module.
selected_map <- bitr(unique(new_selected_hubs$gene_name), fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
universe_symbols <- read.csv(file.path(table_dir, "02_module_assignments_non_grey.csv"), check.names = FALSE)$gene_name
universe_map <- bitr(unique(universe_symbols), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go <- enrichGO(
  gene = unique(selected_map$ENTREZID), universe = unique(universe_map$ENTREZID),
  OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
  pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE
)
go_df <- as.data.frame(go)
write.csv(go_df, file.path(table_dir, "18_selected_module_GO_BP_enrichment.csv"), row.names = FALSE)

gmt_combined <- bind_rows(read.gmt(gmt_legacy_file), read.gmt(gmt_medicus_file)) %>% distinct(term, gene)
kk <- enricher(
  gene = unique(selected_map$ENTREZID), universe = unique(universe_map$ENTREZID),
  TERM2GENE = gmt_combined, pvalueCutoff = 1, qvalueCutoff = 1,
  pAdjustMethod = "BH", minGSSize = 1, maxGSSize = 1000
)
kegg_df <- as.data.frame(kk)
write.csv(kegg_df, file.path(table_dir, "19_selected_module_KEGG_enrichment.csv"), row.names = FALSE)

dot_enrich <- function(df, title, top_n = 5) {
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = paste(title, "(no enriched terms)")))
  d <- df %>% arrange(p.adjust, desc(Count)) %>% slice_head(n = top_n) %>%
    mutate(GeneRatio_num = vapply(strsplit(GeneRatio, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1)),
           Description = factor(Description, levels = rev(Description)))
  ggplot(d, aes(GeneRatio_num, Description, size = Count, color = p.adjust)) +
    geom_point() + scale_color_gradient(low = "#F8766D", high = "#00BFC4") +
    scale_size_continuous(range = c(3, 8)) +
    labs(x = "GeneRatio", y = NULL, title = title) + theme_bw(base_size = 11)
}

p_go <- dot_enrich(go_df, paste0(selected_module, " GO biological processes"), 5)
p_kegg <- dot_enrich(kegg_df, paste0(selected_module, " KEGG pathways"), 5)
ggsave(file.path(figure_dir, "FigureS5E_selected_module_GO_BP.pdf"), p_go, width = 8, height = 4.5, bg = "white")
ggsave(file.path(figure_dir, "FigureS5E_selected_module_GO_BP.png"), p_go, width = 8, height = 4.5, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "FigureS5F_selected_module_KEGG.pdf"), p_kegg, width = 8, height = 4.5, bg = "white")
ggsave(file.path(figure_dir, "FigureS5F_selected_module_KEGG.png"), p_kegg, width = 8, height = 4.5, dpi = 400, bg = "white")

saveRDS(list(
  selected_module = selected_module,
  independent_genes = independent_genes,
  module_correlations = cor_table,
  eligible_sample_data = eligible,
  state_summary = state_summary,
  group_model = group_lm,
  GO = go,
  KEGG = kk
), file.path(output_dir, "downstream_statistics_and_enrichment.rds"))

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_downstream.txt"))
cat("End:", format(Sys.time()), "\n")
