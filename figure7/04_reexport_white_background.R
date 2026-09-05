user_lib <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(patchwork)
  library(ggplot2)
})

output_dir <- "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final/hdWGCNA_rerun_20260824"
figure_dir <- file.path(output_dir, "figures")
seurat_obj <- readRDS(file.path(output_dir, "hdWGCNA_object_final.rds"))

soft_plots <- PlotSoftPowers(seurat_obj)
p_a <- wrap_plots(soft_plots, ncol = 2)
ggsave(file.path(figure_dir, "Figure7A_soft_power.pdf"), p_a, width = 9, height = 4.5, bg = "white")
ggsave(file.path(figure_dir, "Figure7A_soft_power.png"), p_a, width = 9, height = 4.5, dpi = 400, bg = "white")

p_c <- PlotKMEs(seurat_obj, ncol = 3, n_hubs = 20)
ggsave(file.path(figure_dir, "Figure7C_module_hub_genes.pdf"), p_c, width = 12, height = 8, bg = "white")
ggsave(file.path(figure_dir, "Figure7C_module_hub_genes.png"), p_c, width = 12, height = 8, dpi = 400, bg = "white")

feature_plots <- ModuleFeaturePlot(seurat_obj, features = "hMEs", order = TRUE)
p_s5a <- wrap_plots(feature_plots, ncol = 3)
ggsave(file.path(figure_dir, "FigureS5A_module_feature_plots.pdf"), p_s5a, width = 12, height = 8, bg = "white")
ggsave(file.path(figure_dir, "FigureS5A_module_feature_plots.png"), p_s5a, width = 12, height = 8, dpi = 400, bg = "white")

mods <- intersect(c("Fibroblasts1", "Fibroblasts2", "Fibroblasts6"), colnames(seurat_obj[[]]))
p_s5d_list <- lapply(mods, function(mod) {
  VlnPlot(seurat_obj, features = mod, group.by = "Fibro_subtype", pt.size = 0) +
    geom_boxplot(width = 0.25, fill = "white", outlier.shape = NA) +
    xlab("") + ylab("harmonized module eigengene") + NoLegend() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
})
p_s5d <- wrap_plots(p_s5d_list, ncol = 1)
ggsave(file.path(figure_dir, "FigureS5D_modules_1_2_6_by_revised_fibroblast_state.pdf"),
       p_s5d, width = 8, height = 4 * length(p_s5d_list), bg = "white")
ggsave(file.path(figure_dir, "FigureS5D_modules_1_2_6_by_revised_fibroblast_state.png"),
       p_s5d, width = 8, height = 4 * length(p_s5d_list), dpi = 400, bg = "white")

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_white_background_reexport.txt"))
