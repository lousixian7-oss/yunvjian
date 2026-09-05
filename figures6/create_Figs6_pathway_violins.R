#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(patchwork)
})

input_object <- "G:/1Yunvjian/0a26.7.7singlecell/cellchat8.14/cellchat_allcell_original_style_ECM/CellChat_PD_allcell_revised.rds"
output_dir <- "."

cellchat_PD <- readRDS(input_object)

pathway_delta <- c(
  COLLAGEN = 0.603,
  FN1 = 0.160,
  LAMININ = 0.134,
  PERIOSTIN = 0.0513
)

panel_height <- c(
  COLLAGEN = 28,
  FN1 = 14,
  LAMININ = 28,
  PERIOSTIN = 8
)

missing_pathways <- setdiff(names(pathway_delta), cellchat_PD@netP$pathways)
if (length(missing_pathways) > 0) {
  stop("Pathways absent from the PD CellChat object: ", paste(missing_pathways, collapse = ", "))
}

format_delta <- function(pathway) {
  if (pathway == "PERIOSTIN") "0.0513" else sprintf("%.3f", pathway_delta[[pathway]])
}

make_panel <- function(pathway, panel_label) {
  p <- plotGeneExpression(
    cellchat_PD,
    signaling = pathway,
    enriched.only = FALSE
  )

  p <- p &
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.3),
      axis.text.y = element_text(size = 7),
      axis.ticks = element_line(linewidth = 0.25),
      axis.line = element_line(linewidth = 0.35),
      legend.title = element_blank(),
      legend.text = element_text(size = 6.5),
      legend.key.height = grid::unit(0.28, "cm"),
      legend.key.width = grid::unit(0.28, "cm"),
      plot.margin = margin(6, 6, 6, 6)
    )

  p
}

panels <- Map(
  make_panel,
  names(pathway_delta),
  LETTERS[seq_along(pathway_delta)]
)

for (i in seq_along(panels)) {
  pathway <- names(pathway_delta)[[i]]
  ggsave(
    filename = file.path(output_dir, sprintf("Figs6_%s_panel.pdf", pathway)),
    plot = panels[[i]],
    width = 8.0,
    height = panel_height[[pathway]],
    units = "in",
    device = cairo_pdf
  )
}

writeLines(
  c(
    paste0("Input CellChat object: ", input_object),
    paste0("CellChat version: ", as.character(packageVersion("CellChat"))),
    "Displayed condition: PD (matching the original Fig. S6F plotGeneExpression workflow)",
    "Communication-probability differences shown only as title annotations:",
    sprintf("  %s: PD - NC = %s", names(pathway_delta), vapply(names(pathway_delta), format_delta, character(1))),
    "CellChat groups:",
    paste0("  ", levels(cellchat_PD@idents))
  ),
  con = file.path(output_dir, "Figs6_pathway_violins_provenance.txt")
)

message("Created four standalone pathway panel PDFs")
