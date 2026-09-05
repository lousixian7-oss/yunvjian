user_lib <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

output_dir <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final/hdWGCNA_rerun_20260824"
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
ref_dir <- file.path(output_dir, "reference_inputs")
log_con <- file(file.path(output_dir, "03_resume_figures_enrichment.log"), open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Start:", format(Sys.time()), "\n")
cor_table <- read.csv(file.path(table_dir, "06_AGT_module_spearman_all_modules.csv"), check.names = FALSE)
selected_module <- cor_table$module[which.max(cor_table$rho_AGT_mean)]
mods <- cor_table$module
eligible <- read.csv(file.path(table_dir, "05_sample_level_AGT_PI16_module_hMEs.csv"), check.names = FALSE) %>%
  filter(n_fibro >= 50)
state_summary <- read.csv(file.path(table_dir, "07_sample_state_module_hMEs_n10.csv"), check.names = FALSE)
state_global <- read.csv(file.path(table_dir, "08_module_state_kruskal_tests.csv"), check.names = FALSE)
group_lm_coef <- read.csv(file.path(table_dir, "10_selected_module_dataset_adjusted_group_model.csv"), check.names = FALSE)
group_test <- read.csv(file.path(table_dir, "11_selected_module_PD_NC_wilcoxon.csv"), check.names = FALSE)
score_df <- read.csv(file.path(table_dir, "15_independent_module_score_by_sample_state_n10.csv"), check.names = FALSE)
score_kw <- read.csv(file.path(table_dir, "17_independent_score_state_kruskal.csv"), check.names = FALSE)
hub <- read.csv(file.path(table_dir, "03_hub_genes_top100_each_module.csv"), check.names = FALSE)
selected_hubs <- hub %>% filter(module == selected_module) %>% arrange(desc(kME))

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
sig_df <- module_long %>% group_by(module) %>% dplyr::slice(1) %>% ungroup()

p_d <- ggplot(module_long, aes(state_final, module, fill = mean_hME_z)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = sig_df, aes(x = 5.45, y = module, label = sig), inherit.aes = FALSE, hjust = 0, size = 4) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0) +
  coord_cartesian(xlim = c(0.5, 5.8), clip = "off") +
  labs(x = NULL, y = NULL, fill = "mean z-scored hME",
       title = "Module activity across revised fibroblast states",
       subtitle = "GSM x state means (at least 10 cells); stars: BH-adjusted Kruskal-Wallis tests") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(), plot.margin = margin(5.5, 20, 5.5, 5.5))
ggsave(file.path(figure_dir, "Figure7D_sample_level_module_state_heatmap.pdf"), p_d, width = 9, height = 5.5, bg = "white")
ggsave(file.path(figure_dir, "Figure7D_sample_level_module_state_heatmap.png"), p_d, width = 9, height = 5.5, dpi = 400, bg = "white")

sel_cor <- cor_table %>% filter(module == selected_module)
p_e <- ggplot(eligible, aes(x = .data[[selected_module]], y = AGT_mean, color = dataset, shape = group)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "grey35", linewidth = 0.7) +
  labs(x = paste0(selected_module, " mean harmonized eigengene"), y = "Mean AGT expression",
       title = "Sample-level association between AGT and the selected module",
       subtitle = sprintf("Spearman rho = %.2f, p = %.3g, BH-adjusted p = %.3g; GSMs with at least 50 fibroblasts",
                          sel_cor$rho_AGT_mean, sel_cor$p_AGT_mean, sel_cor$p_adjust_AGT_mean)) +
  theme_bw(base_size = 11)
ggsave(file.path(figure_dir, "Figure7E_AGT_selected_module_sample_correlation.pdf"), p_e, width = 7, height = 5, bg = "white")
ggsave(file.path(figure_dir, "Figure7E_AGT_selected_module_sample_correlation.png"), p_e, width = 7, height = 5, dpi = 400, bg = "white")

p_f <- ggplot(score_df, aes(x = factor(state_final, levels = state_levels), y = independent_module_score, fill = state_final)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = dataset), width = 0.14, size = 2, alpha = 0.85) +
  labs(x = NULL, y = "Independent module score", fill = "state", color = "dataset",
       title = paste0(selected_module, " program score across fibroblast states"),
       subtitle = sprintf("Top 25 hubs excluding AGT, PI16 and technical genes; Kruskal-Wallis p = %.3g", score_kw$p)) +
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
       subtitle = sprintf("Sample-level Wilcoxon p = %.3g; dataset-adjusted model p = %.3g\nAt least 50 fibroblasts per GSM",
                          group_test$wilcoxon_p, lm_group_p), color = "dataset", fill = "group") +
  theme_bw(base_size = 11)
ggsave(file.path(figure_dir, "Figure7G_selected_module_PD_NC_sample_level.pdf"), p_g, width = 6.5, height = 5, bg = "white")
ggsave(file.path(figure_dir, "Figure7G_selected_module_PD_NC_sample_level.png"), p_g, width = 6.5, height = 5, dpi = 400, bg = "white")

# Follow the original enrichment code: top 100 hub genes, SYMBOL-to-ENTREZ
# conversion, GO enrichment and combined KEGG legacy/medicus GMT enrichment.
selected_map <- bitr(unique(selected_hubs$gene_name), fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
go <- enrichGO(gene = unique(selected_map$ENTREZID), OrgDb = org.Hs.eg.db,
               pvalueCutoff = 0.05, qvalueCutoff = 1, ont = "all", readable = TRUE)
go_df <- as.data.frame(go)
write.csv(go_df, file.path(table_dir, "18_selected_module_GO_all_enrichment.csv"), row.names = FALSE)

gmt_combined <- bind_rows(
  read.gmt(file.path(ref_dir, "kegg_legacy_2026_1.gmt")),
  read.gmt(file.path(ref_dir, "kegg_medicus_2026_1.gmt"))
) %>% distinct(term, gene)
kk <- enricher(gene = unique(selected_map$ENTREZID), TERM2GENE = gmt_combined,
               pvalueCutoff = 0.5, qvalueCutoff = 1, pAdjustMethod = "fdr",
               minGSSize = 1, maxGSSize = 1000)
kegg_df <- as.data.frame(kk)
write.csv(kegg_df, file.path(table_dir, "19_selected_module_KEGG_enrichment.csv"), row.names = FALSE)

pretty_kegg <- function(x) {
  x <- sub("^KEGG_MEDICUS_REFERENCE_", "", x)
  x <- sub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("ANXA2 S100A10", "ANXA2-S100A10", x)
  x <- gsub("ITGA B", "ITGA/B", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("Anxa2-S100a10", "ANXA2-S100A10", x)
  x <- gsub("Itga/B", "ITGA/B", x)
  x <- gsub("Hcm$", "HCM", x)
  vapply(x, function(s) paste(strwrap(s, width = 42), collapse = "\n"), character(1))
}

dot_enrich <- function(df, title, top_n = 5, bp_only = FALSE, kegg_labels = FALSE) {
  if (bp_only && "ONTOLOGY" %in% names(df)) df <- df %>% filter(ONTOLOGY == "BP")
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = paste(title, "(no enriched terms)")))
  d <- df %>% arrange(p.adjust, desc(Count)) %>% dplyr::slice_head(n = top_n) %>%
    mutate(GeneRatio_num = vapply(strsplit(GeneRatio, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1)),
           plot_label = if (kegg_labels) pretty_kegg(Description) else Description,
           plot_label = factor(plot_label, levels = rev(plot_label)))
  ggplot(d, aes(GeneRatio_num, plot_label, size = Count, color = p.adjust)) +
    geom_point() + scale_color_gradient(low = "#F8766D", high = "#00BFC4") +
    scale_size_continuous(range = c(3, 8)) +
    labs(x = "GeneRatio", y = NULL, title = title) + theme_bw(base_size = 11)
}

p_go <- dot_enrich(go_df, paste0(selected_module, " GO terms"), 5, bp_only = TRUE)
p_kegg <- dot_enrich(kegg_df, paste0(selected_module, " KEGG pathways"), 5, kegg_labels = TRUE)
ggsave(file.path(figure_dir, "FigureS5E_selected_module_GO_BP.pdf"), p_go, width = 8, height = 4.5, bg = "white")
ggsave(file.path(figure_dir, "FigureS5E_selected_module_GO_BP.png"), p_go, width = 8, height = 4.5, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "FigureS5F_selected_module_KEGG.pdf"), p_kegg, width = 9, height = 4.8, bg = "white")
ggsave(file.path(figure_dir, "FigureS5F_selected_module_KEGG.png"), p_kegg, width = 9, height = 4.8, dpi = 400, bg = "white")

writeLines(c(
  paste0("selected_module=", selected_module),
  paste0("AGT_mean_rho=", sel_cor$rho_AGT_mean),
  paste0("AGT_mean_p=", sel_cor$p_AGT_mean),
  paste0("AGT_mean_BH=", sel_cor$p_adjust_AGT_mean),
  paste0("PD_NC_Wilcoxon_p=", group_test$wilcoxon_p),
  paste0("PD_NC_dataset_adjusted_p=", lm_group_p),
  paste0("score_state_KW_p=", score_kw$p)
), file.path(output_dir, "key_results.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_resume_figures_enrichment.txt"))
cat("End:", format(Sys.time()), "\n")
