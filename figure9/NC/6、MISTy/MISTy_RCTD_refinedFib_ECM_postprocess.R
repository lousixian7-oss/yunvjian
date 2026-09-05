#!/usr/bin/env Rscript

.libPaths(c("G:/gurobi/Rlib45", .libPaths()))
suppressPackageStartupMessages({
  library(mistyR)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

source_file <- "MISTy_RCTD_refinedFib_ECM_full.R"
collect_file <- "MISTy_RCTD_refinedFib_ECM_collect.rds"
spatial_file <- "stRNA_RCTD_refinedFib_ECM.rds"
weights_file <- "RCTD_refinedFib_ECM_weights_full.csv"
stopifnot(file.exists(collect_file), file.exists(spatial_file), file.exists(weights_file))

misty_collect <- readRDS(collect_file)

write.csv(misty_collect$improvements, "MISTy_improvements.csv", row.names = FALSE)
write.csv(misty_collect$contributions, "MISTy_contributions.csv", row.names = FALSE)
write.csv(misty_collect$importances.aggregated, "MISTy_importances_aggregated.csv", row.names = FALSE)

save_misty_plot <- function(filename, expression, width, height) {
  grDevices::pdf(filename, width = width, height = height, onefile = TRUE)
  tryCatch(force(expression), finally = grDevices::dev.off())
  invisible(NULL)
}

save_misty_plot(
  "MISTy_RCTD_refinedFib_ECM_improvement.pdf",
  plot_improvement_stats(misty_collect),
  width = 9,
  height = 6
)

view_names <- unique(misty_collect$importances.aggregated$view)
for (view_name in view_names) {
  safe_view <- gsub("[^A-Za-z0-9_.-]", "_", view_name)
  save_misty_plot(
    paste0("MISTy_RCTD_refinedFib_ECM_interactions_", safe_view, ".pdf"),
    plot_interaction_heatmap(misty_collect, view = view_name, cutoff = 0),
    width = 9,
    height = 8
  )
}

stRNA <- readRDS(spatial_file)
weights_df <- read.csv(weights_file, check.names = FALSE, stringsAsFactors = FALSE)
rownames(weights_df) <- weights_df$spot_id
weights_df$spot_id <- NULL
rctd_mat <- as.matrix(weights_df)
rctd_mat <- rctd_mat / rowSums(rctd_mat)
colnames(rctd_mat) <- make.names(colnames(rctd_mat), unique = TRUE)

images <- attr(stRNA, "images", exact = TRUE)
image_attrs <- attributes(images[[1]])
centroid_attrs <- attributes(image_attrs$boundaries[["centroids"]])
coords <- as.data.frame(centroid_attrs$coords)
rownames(coords) <- centroid_attrs$cells
colnames(coords)[1:2] <- c("x", "y")
coords <- coords[, c("x", "y"), drop = FALSE]

common_spots <- intersect(rownames(rctd_mat), rownames(coords))
plot_df <- cbind(
  coords[common_spots, , drop = FALSE],
  as.data.frame(rctd_mat[common_spots, , drop = FALSE])
)
focus_long <- plot_df %>%
  mutate(spot_id = rownames(plot_df)) %>%
  select(spot_id, x, y, ECM_Fib, Activated_Fib) %>%
  pivot_longer(c(ECM_Fib, Activated_Fib), names_to = "celltype", values_to = "proportion")

p_focus <- ggplot(focus_long, aes(x = x, y = y, color = proportion)) +
  geom_point(size = 1.2, alpha = 0.95) +
  scale_color_viridis_c(option = "turbo") +
  scale_y_reverse() +
  coord_fixed() +
  facet_wrap(~celltype, nrow = 1) +
  theme_void(base_size = 12) +
  theme(legend.position = "top", strip.text = element_text(face = "bold")) +
  labs(color = "RCTD proportion")
ggsave("ECM_Fib_Activated_Fib_spatial_MISTy_input.pdf", p_focus, width = 10, height = 5)

cat("MISTy postprocessing completed. Views:", paste(view_names, collapse = ", "), "\n")
