user_lib <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(igraph)
})

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 7)

input_rds <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final/results/objects/fibroblast_revision_final_no_unassigned.rds"
output_dir <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final/hdWGCNA_rerun_20260824"
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "TOM"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "ModuleNetworks"), recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

log_con <- file(file.path(output_dir, "01_source_faithful_hdwgcna.log"), open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Start:", format(Sys.time()), "\n")
cat("Input:", input_rds, "\n")
cat("Output:", output_dir, "\n")
cat("R:", as.character(getRversion()), "\n")
for (pkg in c("Seurat", "SeuratObject", "hdWGCNA", "WGCNA", "UCell")) {
  cat(pkg, as.character(packageVersion(pkg)), "\n")
}

yjsl <- readRDS(input_rds)
required_meta <- c("orig.ident", "state_final", "cellType", "dataset", "group")
stopifnot(all(required_meta %in% colnames(yjsl[[]])))
stopifnot("harmony" %in% Reductions(yjsl))

# Source substitution only: the revised fibroblast state is exposed under the
# original script's Fibro_subtype field, while all original hdWGCNA parameters
# and the original harmony reduction are retained.
yjsl$Fibro_subtype <- factor(
  as.character(yjsl$state_final),
  levels = c("PI16_Fib", "Adventitial_Fib", "ECM_Fib", "Activated_Fib", "Inflammatory_Fib")
)
yjsl$cellType <- factor("Fibroblasts", levels = "Fibroblasts")
Idents(yjsl) <- yjsl$Fibro_subtype
DefaultAssay(yjsl) <- "RNA"

write.csv(
  as.data.frame(table(
    orig.ident = yjsl$orig.ident,
    Fibro_subtype = yjsl$Fibro_subtype,
    dataset = yjsl$dataset,
    group = yjsl$group
  )),
  file.path(table_dir, "00_input_cells_by_sample_state_dataset_group.csv"),
  row.names = FALSE
)

seurat_obj <- yjsl
rm(yjsl)
gc()

# This is the executable correction of the original line-order typo: the
# original file removed yjsl and then attempted to subset yjsl. The intended
# object is seurat_obj; all filtering semantics are otherwise unchanged.
fibro_wgcna <- subset(seurat_obj, subset = !is.na(Fibro_subtype))
rm(seurat_obj)
gc()
print(table(fibro_wgcna$Fibro_subtype))

fibro_wgcna <- SetupForWGCNA(
  fibro_wgcna,
  gene_select = "fraction",
  fraction = 0.03,
  wgcna_name = "Fibroblast"
)

fibro_wgcna <- MetacellsByGroups(
  seurat_obj = fibro_wgcna,
  group.by = c("orig.ident", "Fibro_subtype"),
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  ident.group = "Fibro_subtype"
)

fibro_wgcna <- NormalizeMetacells(fibro_wgcna)

fibro_wgcna <- SetDatExpr(
  fibro_wgcna,
  group_name = NULL,
  group.by = NULL,
  assay = "RNA",
  slot = "data"
)

fibro_wgcna <- TestSoftPowers(fibro_wgcna, networkType = "signed")
seurat_obj <- fibro_wgcna
rm(fibro_wgcna)
gc()

power_table <- GetPowerTable(seurat_obj)
write.csv(power_table, file.path(table_dir, "01_soft_power_table.csv"), row.names = FALSE)
plot_list <- PlotSoftPowers(seurat_obj)
p1 <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(figure_dir, "Figure7A_soft_power.pdf"), p1, width = 9, height = 4.5)
ggsave(file.path(figure_dir, "Figure7A_soft_power.png"), p1, width = 9, height = 4.5, dpi = 400)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 5,
  setDatExpr = FALSE,
  tom_name = "Fibroblast",
  overwrite_tom = TRUE
)

# Persist the completed network before any downstream plotting. This checkpoint
# prevents a graphics-device failure from discarding the expensive TOM/network.
saveRDS(seurat_obj, file.path(output_dir, "hdWGCNA_constructed_network_checkpoint.rds"), compress = FALSE)

# PlotDendrogram draws directly to the active graphics device and returns a
# numeric value, so it must not be passed to ggsave().
pdf(file.path(figure_dir, "Figure7B_dendrogram.pdf"), width = 10, height = 5)
PlotDendrogram(seurat_obj, main = "Fibroblasts hdWGCNA Dendrogram")
dev.off()
png(file.path(figure_dir, "Figure7B_dendrogram.png"), width = 4000, height = 2000, res = 400)
PlotDendrogram(seurat_obj, main = "Fibroblasts hdWGCNA Dendrogram")
dev.off()

seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "orig.ident")

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "cellType",
  group_name = "Fibroblasts"
)

seurat_obj <- ResetModuleNames(seurat_obj, new_name = "Fibroblasts")

p3 <- PlotKMEs(seurat_obj, ncol = 3, n_hubs = 20)
ggsave(file.path(figure_dir, "Figure7C_module_hub_genes.pdf"), p3, width = 12, height = 8)
ggsave(file.path(figure_dir, "Figure7C_module_hub_genes.png"), p3, width = 12, height = 8, dpi = 400)

modules <- GetModules(seurat_obj) %>% subset(module != "grey")
hub_df <- GetHubGenes(seurat_obj, n_hubs = 100)
write.csv(modules, file.path(table_dir, "02_module_assignments_non_grey.csv"), row.names = FALSE)
write.csv(hub_df, file.path(table_dir, "03_hub_genes_top100_each_module.csv"), row.names = FALSE)
write.csv(as.data.frame(table(GetModules(seurat_obj)$module)),
          file.path(table_dir, "04_module_sizes.csv"), row.names = FALSE)

saveRDS(seurat_obj, file.path(output_dir, "hdWGCNA_network_checkpoint.rds"), compress = FALSE)

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method = "Seurat"
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = "hMEs",
  order = TRUE
)
p4 <- wrap_plots(plot_list, ncol = 3)
ggsave(file.path(figure_dir, "FigureS5A_module_feature_plots.pdf"), p4, width = 12, height = 8)
ggsave(file.path(figure_dir, "FigureS5A_module_feature_plots.png"), p4, width = 12, height = 8, dpi = 400)

# ModuleCorrelogram also draws to the active graphics device in hdWGCNA 0.4.12.
pdf(file.path(figure_dir, "FigureS5B_module_correlogram.pdf"), width = 8, height = 7)
ModuleCorrelogram(seurat_obj)
dev.off()
png(file.path(figure_dir, "FigureS5B_module_correlogram.png"), width = 3200, height = 2800, res = 400)
ModuleCorrelogram(seurat_obj)
dev.off()

MEs <- GetMEs(seurat_obj, harmonized = TRUE)
mods <- setdiff(colnames(MEs), "grey")
duplicate_me_cols <- intersect(colnames(seurat_obj[[]]), colnames(MEs))
if (length(duplicate_me_cols) > 0) {
  seurat_obj@meta.data[, duplicate_me_cols] <- NULL
}
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p6 <- DotPlot(seurat_obj, features = mods, group.by = "Fibro_subtype") +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
ggsave(file.path(figure_dir, "Figure7D_modules_by_revised_fibroblast_state.pdf"), p6, width = 9, height = 6)
ggsave(file.path(figure_dir, "Figure7D_modules_by_revised_fibroblast_state.png"), p6, width = 9, height = 6, dpi = 400)

vln_mods <- intersect(c("Fibroblasts1", "Fibroblasts2", "Fibroblasts6"), mods)
if (length(vln_mods) > 0) {
  p9_list <- lapply(vln_mods, function(mod) {
    VlnPlot(seurat_obj, features = mod, group.by = "Fibro_subtype", pt.size = 0) +
      geom_boxplot(width = 0.25, fill = "white", outlier.shape = NA) +
      xlab("") + ylab("harmonized module eigengene") + NoLegend() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })
  p9 <- wrap_plots(p9_list, ncol = 1)
  ggsave(file.path(figure_dir, "FigureS5D_modules_1_2_6_by_revised_fibroblast_state.pdf"),
         p9, width = 8, height = 4 * length(p9_list))
  ggsave(file.path(figure_dir, "FigureS5D_modules_1_2_6_by_revised_fibroblast_state.png"),
         p9, width = 8, height = 4 * length(p9_list), dpi = 400)
}

pdf(file.path(figure_dir, "FigureS5C_hub_gene_network.pdf"), width = 9, height = 6)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2,
  n_other = 5,
  edge_prop = 0.55,
  mods = "all"
)
dev.off()

saveRDS(seurat_obj, file.path(output_dir, "hdWGCNA_object_final.rds"), compress = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_hdwgcna_rerun.txt"))
cat("End:", format(Sys.time()), "\n")
