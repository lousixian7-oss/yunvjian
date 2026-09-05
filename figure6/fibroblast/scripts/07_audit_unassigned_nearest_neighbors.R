args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1]] else "G:/1Yunvjian/0a26.7.7singlecell/8.10/fibroblast_revision_final"

suppressPackageStartupMessages({
  library(Seurat)
  library(RANN)
})

obj <- readRDS(file.path(out_dir, "results/objects/fibroblast_revision_final.rds"))
emb <- Embeddings(obj, "harmony_fibro_final")
lab <- as.character(obj$state_final)
query_idx <- which(lab == "Unassigned_Fib")
ref_idx <- which(lab != "Unassigned_Fib")

nn <- RANN::nn2(
  data = emb[ref_idx, , drop = FALSE],
  query = emb[query_idx, , drop = FALSE],
  k = 30
)$nn.idx
neighbor_labels <- matrix(lab[ref_idx][nn], nrow = nrow(nn), ncol = ncol(nn))
vote_counts <- sort(table(as.vector(neighbor_labels)), decreasing = TRUE)

per_cell_majority <- apply(neighbor_labels, 1, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})
majority_counts <- sort(table(per_cell_majority), decreasing = TRUE)

centroids <- rowsum(emb[ref_idx, , drop = FALSE], lab[ref_idx]) /
  as.vector(table(lab[ref_idx])[rownames(rowsum(emb[ref_idx, , drop = FALSE], lab[ref_idx]))])
query_centroid <- colMeans(emb[query_idx, , drop = FALSE])
centroid_distance <- sort(apply(centroids, 1, function(x) sqrt(sum((x - query_centroid)^2))))

audit <- data.frame(
  metric = c(
    paste0("30NN_all_votes_", names(vote_counts)),
    paste0("30NN_cell_majority_", names(majority_counts)),
    paste0("centroid_distance_", names(centroid_distance))
  ),
  value = c(as.numeric(vote_counts), as.numeric(majority_counts), as.numeric(centroid_distance))
)

write.csv(
  audit,
  file.path(out_dir, "results/tables/figure6_revised_GK/unassigned_assignment_audit.csv"),
  row.names = FALSE
)
print(audit)
