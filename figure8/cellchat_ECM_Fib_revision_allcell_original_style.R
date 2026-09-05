library(Seurat)
library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(future)


user_lib <- "C:/Users/32266/AppData/Local/R/win-library/4.5"

if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}


out_dir <- "cellchat_allcell_original_style_ECM"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


set.seed(20260822)


cell_keep <- c(
  "T cells",
  "B cells",
  "Macrophages",
  "Endothelial cells",
  "Plasma cells",
  "Epithelial cells",
  "Pericytes",
  "Mast cells",
  "Neutrophils",
  "pDCs",
  "PI16_Fib",
  "Activated_Fib",
  "ECM_Fib",
  "Inflammatory_Fib",
  "Adventitial_Fib"
)


NC_color <- "#00BFC4"
PD_color <- "#F8766D"


saved_corrected_dir <-
  "cellchat_revised_fib_state_final"


saved_NC <- file.path(
  saved_corrected_dir,
  "CellChat_NC_revised_fib_matched.rds"
)


saved_PD <- file.path(
  saved_corrected_dir,
  "CellChat_PD_revised_fib_matched.rds"
)


stopifnot(
  file.exists(saved_NC),
  file.exists(saved_PD)
)


cellchat_NC <- readRDS(
  saved_NC
)


cellchat_PD <- readRDS(
  saved_PD
)


stopifnot(
  identical(
    rownames(cellchat_NC@net$weight),
    rownames(cellchat_PD@net$weight)
  )
)


cellchat_merge <- mergeCellChat(
  object.list = list(
    NC = cellchat_NC,
    PD = cellchat_PD
  ),
  add.names = c(
    "NC",
    "PD"
  )
)


saveRDS(
  cellchat_NC,
  file.path(
    out_dir,
    "cellchat_NC_allcell_revised.rds"
  )
)


saveRDS(
  cellchat_PD,
  file.path(
    out_dir,
    "CellChat_PD_allcell_revised.rds"
  )
)


saveRDS(
  cellchat_merge,
  file.path(
    out_dir,
    "cellchat_merge_NC_PD_allcell_revised.rds"
  )
)


all_types <- rownames(
  cellchat_NC@net$weight
)


celltype_colors <- setNames(
  scales::hue_pal()(length(all_types)),
  all_types
)


celltype_colors["ECM_Fib"] <- "#2A9D8F"
celltype_colors["Activated_Fib"] <- "#E64B35"


########## Figure 8A: 通讯数量和强度比较 ##########


p_count <- compareInteractions(
  cellchat_merge,
  show.legend = FALSE,
  group = c(1, 2),
  color.use = c(
    NC_color,
    PD_color
  )
) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )


p_strength <- compareInteractions(
  cellchat_merge,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight",
  color.use = c(
    NC_color,
    PD_color
  )
) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )


p_overall <- p_count + p_strength +
  plot_layout(ncol = 2)


ggsave(
  file.path(
    out_dir,
    "Fig8A_allcell_communication_number_strength.pdf"
  ),
  p_overall,
  width = 6,
  height = 4
)


ggsave(
  file.path(
    out_dir,
    "Fig8A_allcell_communication_number_strength.png"
  ),
  p_overall,
  width = 6,
  height = 4,
  dpi = 600
)


########## Figure 8B: 全部细胞 incoming 和 outgoing ##########


get_strength <- function(cellchat_obj, all_types) {

  weight_mat <- cellchat_obj@net$weight

  outgoing <- rowSums(weight_mat)
  incoming <- colSums(weight_mat)

  data.frame(
    CellType = all_types,
    Incoming = incoming[all_types],
    Outgoing = outgoing[all_types]
  )
}


strength_NC <- get_strength(
  cellchat_NC,
  all_types
)


strength_PD <- get_strength(
  cellchat_PD,
  all_types
)


write.csv(
  bind_rows(
    mutate(strength_NC, Condition = "NC"),
    mutate(strength_PD, Condition = "PD")
  ),
  file.path(
    out_dir,
    "allcell_incoming_outgoing_strength.csv"
  ),
  row.names = FALSE
)


p_role_NC <- ggplot(
  strength_NC,
  aes(
    x = Incoming,
    y = Outgoing,
    color = CellType,
    fill = CellType
  )
) +
  geom_point(
    size = 5,
    shape = 21,
    stroke = 0.8
  ) +
  geom_text_repel(
    aes(label = CellType),
    size = 3.2,
    max.overlaps = Inf,
    seed = 20260822
  ) +
  scale_color_manual(
    values = celltype_colors
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  labs(
    title = "NC",
    x = "Incoming interaction strength",
    y = "Outgoing interaction strength"
  )


p_role_PD <- ggplot(
  strength_PD,
  aes(
    x = Incoming,
    y = Outgoing,
    color = CellType,
    fill = CellType
  )
) +
  geom_point(
    size = 5,
    shape = 21,
    stroke = 0.8
  ) +
  geom_text_repel(
    aes(label = CellType),
    size = 3.2,
    max.overlaps = Inf,
    seed = 20260822
  ) +
  scale_color_manual(
    values = celltype_colors
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  labs(
    title = "PD",
    x = "Incoming interaction strength",
    y = "Outgoing interaction strength"
  )


x_lim <- range(
  c(
    strength_NC$Incoming,
    strength_PD$Incoming
  )
)


y_lim <- range(
  c(
    strength_NC$Outgoing,
    strength_PD$Outgoing
  )
)


p_role_NC <- p_role_NC +
  coord_cartesian(
    xlim = x_lim,
    ylim = y_lim
  )


p_role_PD <- p_role_PD +
  coord_cartesian(
    xlim = x_lim,
    ylim = y_lim
  )


p_role <- p_role_NC + p_role_PD


ggsave(
  file.path(
    out_dir,
    "Fig8B_allcell_incoming_outgoing_NC_PD.pdf"
  ),
  p_role,
  width = 11,
  height = 5.5
)


ggsave(
  file.path(
    out_dir,
    "Fig8B_allcell_incoming_outgoing_NC_PD.png"
  ),
  p_role,
  width = 11,
  height = 5.5,
  dpi = 600
)


########## Figure 8C: 全部细胞 PD-NC 差异环形图 ##########


weight_diff <- cellchat_PD@net$weight -
  cellchat_NC@net$weight


diff_weight_up <- weight_diff

diff_weight_up[
  diff_weight_up < 0
] <- 0


diff_weight_down <- -weight_diff

diff_weight_down[
  diff_weight_down < 0
] <- 0


group_size <- as.numeric(
  table(cellchat_PD@idents)[all_types]
)


pdf(
  file.path(
    out_dir,
    "Fig8C_PD_enhanced_circle.pdf"
  ),
  width = 7,
  height = 7,
  useDingbats = FALSE
)


netVisual_circle(
  diff_weight_up,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "PD-enhanced interaction strength"
)


dev.off()


pdf(
  file.path(
    out_dir,
    "Fig8C_NC_enhanced_circle.pdf"
  ),
  width = 7,
  height = 7,
  useDingbats = FALSE
)


netVisual_circle(
  diff_weight_down,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  color.use = celltype_colors,
  title.name = "NC-enhanced interaction strength"
)


dev.off()


########## 筛选 ECM_Fib -> Activated_Fib 最可能通路 ##########


ecm_PD <- subsetCommunication(
  cellchat_PD,
  sources.use = "ECM_Fib",
  targets.use = "Activated_Fib"
)


ecm_NC <- subsetCommunication(
  cellchat_NC,
  sources.use = "ECM_Fib",
  targets.use = "Activated_Fib"
)


path_PD <- ecm_PD %>%
  group_by(pathway_name) %>%
  summarise(
    prob_PD = sum(prob),
    n_LR_PD = n(),
    .groups = "drop"
  )


path_NC <- ecm_NC %>%
  group_by(pathway_name) %>%
  summarise(
    prob_NC = sum(prob),
    n_LR_NC = n(),
    .groups = "drop"
  )


pathway_ranking <- full_join(
  path_PD,
  path_NC,
  by = "pathway_name"
) %>%
  mutate(
    across(
      where(is.numeric),
      ~replace_na(.x, 0)
    ),
    diff = prob_PD - prob_NC,
    fold_change =
      (prob_PD + 1e-8) /
      (prob_NC + 1e-8)
  ) %>%
  arrange(
    desc(diff)
  )


write.csv(
  pathway_ranking,
  file.path(
    out_dir,
    "ECM_Fib_to_Activated_Fib_pathway_ranking.csv"
  ),
  row.names = FALSE
)


ecm_priority <- c(
  "COLLAGEN",
  "FN1",
  "LAMININ",
  "PERIOSTIN",
  "THBS",
  "PDGF"
)


selected_pathways <- pathway_ranking %>%
  filter(
    pathway_name %in% ecm_priority,
    diff > 0
  ) %>%
  slice_head(n = 4) %>%
  pull(pathway_name)


stopifnot(
  length(selected_pathways) >= 1
)


writeLines(
  selected_pathways,
  file.path(
    out_dir,
    "selected_ECM_Fib_to_Activated_Fib_pathways.txt"
  )
)


########## Figure 8D: 重点通路全部细胞网络 ##########


pdf(
  file.path(
    out_dir,
    "Fig8D_allcell_selected_pathway_networks.pdf"
  ),
  width = 12,
  height = 12,
  useDingbats = FALSE
)


par(
  mfrow = c(2, 2),
  xpd = TRUE
)


for (pathway_now in selected_pathways) {

  netVisual_aggregate(
    cellchat_PD,
    signaling = pathway_now,
    signaling.name = paste0(
      pathway_now,
      " signaling - PD"
    ),
    layout = "circle",
    color.use = celltype_colors,
    vertex.label.cex = 0.8
  )
}


dev.off()


########## Figure 8E: ECM_Fib -> Activated_Fib 受体配体比较 ##########


ecm_PD_selected <- ecm_PD %>%
  filter(
    pathway_name %in% selected_pathways
  )


ecm_NC_selected <- ecm_NC %>%
  filter(
    pathway_name %in% selected_pathways
  )


compare_df <- full_join(
  ecm_PD_selected,
  ecm_NC_selected,
  by = c(
    "interaction_name",
    "pathway_name",
    "ligand",
    "receptor"
  ),
  suffix = c(
    "_PD",
    "_NC"
  )
) %>%
  mutate(
    prob_PD = replace_na(prob_PD, 0),
    prob_NC = replace_na(prob_NC, 0),
    diff = prob_PD - prob_NC,
    pathway_name = factor(
      pathway_name,
      levels = selected_pathways
    ),
    cell_pair = "ECM_Fib -> Activated_Fib"
  ) %>%
  arrange(
    pathway_name,
    desc(diff)
  )


write.csv(
  compare_df,
  file.path(
    out_dir,
    "ECM_Fib_to_Activated_Fib_selected_LR_PD_vs_NC.csv"
  ),
  row.names = FALSE
)


top_compare_df <- compare_df %>%
  group_by(pathway_name) %>%
  slice_max(
    order_by = abs(diff),
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup()


top_LR <- data.frame(
  interaction_name = unique(
    top_compare_df$interaction_name
  )
)


p_lr_compare <- netVisual_bubble(
  cellchat_merge,
  sources.use = "ECM_Fib",
  targets.use = "Activated_Fib",
  comparison = c(1, 2),
  pairLR.use = top_LR,
  remove.isolate = TRUE,
  title.name =
    "ECM_Fib -> Activated_Fib key ligand-receptor interactions"
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_text(
      size = 8
    ),
    legend.position = "right",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  labs(
    x = NULL,
    y = "Ligand-Receptor pair"
  )


p_lr_diff <- ggplot(
  top_compare_df,
  aes(
    x = interaction_name,
    y = cell_pair,
    size = abs(diff),
    color = diff
  )
) +
  geom_point() +
  facet_wrap(
    ~pathway_name,
    scales = "free_x",
    nrow = 1
  ) +
  scale_color_gradient2(
    low = NC_color,
    mid = "white",
    high = PD_color,
    midpoint = 0
  ) +
  scale_size_continuous(
    range = c(1.5, 7)
  ) +
  theme_classic() +
  theme(
    strip.position = "top",
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      size = 7
    ),
    axis.text.y = element_text(
      size = 9
    ),
    panel.spacing.x = grid::unit(
      0.3,
      "lines"
    ),
    legend.position = "right"
  ) +
  labs(
    x = "Ligand-Receptor pair",
    y = NULL,
    color = "PD - NC",
    size = "|difference|"
  )


p_E <- p_lr_compare / p_lr_diff +
  plot_layout(
    heights = c(1.4, 1)
  )


ggsave(
  file.path(
    out_dir,
    "Fig8E_ECM_Fib_to_Activated_Fib_LR_comparison.pdf"
  ),
  p_E,
  width = 13,
  height = 9
)


ggsave(
  file.path(
    out_dir,
    "Fig8E_ECM_Fib_to_Activated_Fib_LR_comparison.png"
  ),
  p_E,
  width = 13,
  height = 9,
  dpi = 600
)


global_summary <- data.frame(
  Condition = c(
    "NC",
    "PD"
  ),
  Interaction_number = c(
    sum(cellchat_NC@net$count),
    sum(cellchat_PD@net$count)
  ),
  Interaction_strength = c(
    sum(cellchat_NC@net$weight),
    sum(cellchat_PD@net$weight)
  )
)


write.csv(
  global_summary,
  file.path(
    out_dir,
    "CellChat_global_summary_NC_PD.csv"
  ),
  row.names = FALSE
)


writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    out_dir,
    "sessionInfo.txt"
  )
)


message(
  "Completed: ",
  normalizePath(out_dir)
)
