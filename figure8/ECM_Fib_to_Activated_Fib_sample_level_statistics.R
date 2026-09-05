## Sample-level differential statistics for ECM_Fib -> Activated_Fib communication.
##
## This is a biological-replicate sensitivity analysis for Figure 8F. It uses
## the same composition-matched, population.size = FALSE CellChat objects as
## Figure 8, but treats each sample (not each cell) as the replicate.

user_lib <- "C:/Users/32266/AppData/Local/R/win-library/4.5"
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

base_dir <- "G:/1Yunvjian/0a26.7.7singlecell/cellchat8.14"
object_dir <- file.path(base_dir, "cellchat_allcell_original_style_ECM")
out_dir <- file.path(object_dir, "ECM_to_Activated_sample_level_statistics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nc_file <- file.path(object_dir, "cellchat_NC_allcell_revised.rds")
pd_file <- file.path(object_dir, "CellChat_PD_allcell_revised.rds")
stopifnot(file.exists(nc_file), file.exists(pd_file))

sender <- "ECM_Fib"
receiver <- "Activated_Fib"
pathways_tested <- c("COLLAGEN", "FN1", "LAMININ", "PERIOSTIN")
min_cells <- 10L
kh <- 0.5
hill_n <- 1
bootstrap_repetitions <- 10000L
base_seed <- 20260904L

message("Loading composition-matched CellChat objects...")
cellchat_nc <- readRDS(nc_file)
cellchat_pd <- readRDS(pd_file)

stopifnot(
  identical(rownames(cellchat_nc@data), rownames(cellchat_pd@data)),
  all(c("sample", "dataset", "CellChat_annotation") %in% colnames(cellchat_nc@meta)),
  all(c("sample", "dataset", "CellChat_annotation") %in% colnames(cellchat_pd@meta))
)

prepare_one <- function(object, condition) {
  meta <- object@meta[colnames(object@data), , drop = FALSE]
  meta$condition_for_test <- condition
  list(data = object@data, meta = meta)
}

nc <- prepare_one(cellchat_nc, "NC")
pd <- prepare_one(cellchat_pd, "PD")
expression_matrix <- cbind(nc$data, pd$data)
meta <- rbind(nc$meta, pd$meta)
stopifnot(identical(colnames(expression_matrix), rownames(meta)))

sample_map <- meta %>%
  tibble::rownames_to_column("barcode") %>%
  transmute(
    barcode,
    sample = as.character(sample),
    condition = as.character(condition_for_test),
    dataset = as.character(dataset),
    cell_type = as.character(CellChat_annotation)
  )

sample_counts <- sample_map %>%
  filter(cell_type %in% c(sender, receiver)) %>%
  count(sample, condition, dataset, cell_type, name = "n_cells") %>%
  complete(
    nesting(sample, condition, dataset),
    cell_type = c(sender, receiver),
    fill = list(n_cells = 0L)
  ) %>%
  pivot_wider(names_from = cell_type, values_from = n_cells) %>%
  arrange(condition, sample)

write.csv(
  sample_counts,
  file.path(out_dir, "01_sample_cell_counts.csv"),
  row.names = FALSE
)

eligible <- sample_counts %>%
  filter(.data[[sender]] >= min_cells, .data[[receiver]] >= min_cells)

if (n_distinct(eligible$condition) != 2L) {
  stop("Both NC and PD must contain eligible samples.")
}
if (n_distinct(eligible$dataset) != 1L) {
  stop(
    "Eligible samples span multiple datasets. Use a dataset-stratified test instead of the current exact test."
  )
}

write.csv(
  eligible,
  file.path(out_dir, "02_eligible_samples.csv"),
  row.names = FALSE
)

## Restrict the sample-level test to interactions that were inferred for the
## same sender-receiver pair in at least one of the two aggregate objects. This
## aligns the statistical sensitivity analysis with the interactions displayed
## in Figure 8E-F and avoids testing unrelated database entries.
detected_nc <- subsetCommunication(
  cellchat_nc, sources.use = sender, targets.use = receiver
)
detected_pd <- subsetCommunication(
  cellchat_pd, sources.use = sender, targets.use = receiver
)
detected_names <- union(
  detected_nc$interaction_name[detected_nc$pathway_name %in% pathways_tested],
  detected_pd$interaction_name[detected_pd$pathway_name %in% pathways_tested]
)

db_interaction <- CellChatDB.human[["interaction"]]
db_complex <- CellChatDB.human[["complex"]]
lr_table <- db_interaction %>%
  tibble::rownames_to_column("database_row") %>%
  filter(
    interaction_name %in% detected_names,
    pathway_name %in% pathways_tested
  ) %>%
  distinct(interaction_name, pathway_name, ligand, receptor, .keep_all = TRUE) %>%
  arrange(match(pathway_name, pathways_tested), interaction_name)

if (nrow(lr_table) == 0L) stop("No displayed ECM pathway interactions were found in CellChatDB.human.")

complex_subunits <- function(entity) {
  if (!(entity %in% rownames(db_complex))) return(character())
  cols <- grep("^subunit", colnames(db_complex), value = TRUE)
  ans <- unlist(db_complex[entity, cols, drop = FALSE], use.names = FALSE)
  unique(ans[!is.na(ans) & nzchar(ans)])
}

genes_needed <- unique(c(lr_table$ligand, lr_table$receptor))
complex_entities <- intersect(genes_needed, rownames(db_complex))
genes_needed <- unique(c(
  setdiff(genes_needed, complex_entities),
  unlist(lapply(complex_entities, complex_subunits), use.names = FALSE)
))
genes_needed <- intersect(genes_needed, rownames(expression_matrix))
if (length(genes_needed) == 0L) stop("None of the required signaling genes are present in the CellChat data slot.")
expression_matrix <- expression_matrix[genes_needed, , drop = FALSE]

entity_expression <- function(entity, average_expression) {
  if (entity %in% names(average_expression)) {
    return(unname(average_expression[[entity]]))
  }
  subunits <- complex_subunits(entity)
  if (length(subunits) == 0L || !all(subunits %in% names(average_expression))) {
    return(0)
  }
  values <- unname(average_expression[subunits])
  if (any(!is.finite(values)) || any(values <= 0)) return(0)
  exp(mean(log(values)))
}

trimean_by_gene <- function(cells) {
  if (length(cells) == 0L) return(setNames(rep(0, nrow(expression_matrix)), rownames(expression_matrix)))
  values <- as.matrix(expression_matrix[, cells, drop = FALSE])
  ans <- apply(values, 1, CellChat:::triMean)
  ans[!is.finite(ans)] <- 0
  ans
}

message("Computing sample-level CellChat-compatible core probabilities...")
sample_scores <- lapply(seq_len(nrow(eligible)), function(i) {
  sample_i <- eligible$sample[[i]]
  condition_i <- eligible$condition[[i]]
  dataset_i <- eligible$dataset[[i]]
  sender_cells <- sample_map$barcode[
    sample_map$sample == sample_i & sample_map$cell_type == sender
  ]
  receiver_cells <- sample_map$barcode[
    sample_map$sample == sample_i & sample_map$cell_type == receiver
  ]
  sender_avg <- trimean_by_gene(sender_cells)
  receiver_avg <- trimean_by_gene(receiver_cells)

  bind_rows(lapply(seq_len(nrow(lr_table)), function(j) {
    ligand_expression <- entity_expression(lr_table$ligand[[j]], sender_avg)
    receptor_expression <- entity_expression(lr_table$receptor[[j]], receiver_avg)
    mass_action <- ligand_expression * receptor_expression
    probability_core <- mass_action^hill_n / (kh^hill_n + mass_action^hill_n)
    data.frame(
      sample = sample_i,
      condition = condition_i,
      dataset = dataset_i,
      n_sender = eligible[[sender]][[i]],
      n_receiver = eligible[[receiver]][[i]],
      interaction_name = lr_table$interaction_name[[j]],
      pathway_name = lr_table$pathway_name[[j]],
      ligand = lr_table$ligand[[j]],
      receptor = lr_table$receptor[[j]],
      ligand_trimean = ligand_expression,
      receptor_trimean = receptor_expression,
      mass_action = mass_action,
      probability_core = probability_core,
      stringsAsFactors = FALSE
    )
  }))
}) %>% bind_rows()

write.csv(
  sample_scores,
  file.path(out_dir, "03_sample_level_LR_core_probabilities.csv"),
  row.names = FALSE
)

pathway_scores <- sample_scores %>%
  group_by(sample, condition, dataset, n_sender, n_receiver, pathway_name) %>%
  summarise(
    pathway_score = sum(probability_core),
    n_LR = n(),
    .groups = "drop"
  ) %>%
  mutate(pathway_name = factor(pathway_name, levels = pathways_tested)) %>%
  arrange(pathway_name, condition, sample)

write.csv(
  pathway_scores,
  file.path(out_dir, "04_sample_level_pathway_scores.csv"),
  row.names = FALSE
)

exact_permutation_test <- function(values, condition) {
  condition <- factor(condition, levels = c("NC", "PD"))
  n_total <- length(values)
  n_pd <- sum(condition == "PD")
  observed <- mean(values[condition == "PD"]) - mean(values[condition == "NC"])
  combinations <- combn(n_total, n_pd)
  total_sum <- sum(values)
  permuted <- apply(combinations, 2, function(pd_index) {
    pd_sum <- sum(values[pd_index])
    pd_sum / n_pd - (total_sum - pd_sum) / (n_total - n_pd)
  })
  tolerance <- sqrt(.Machine$double.eps)
  list(
    delta_mean = observed,
    p_two_sided = mean(abs(permuted) >= abs(observed) - tolerance),
    p_one_sided_PD_higher = mean(permuted >= observed - tolerance),
    n_permutations = length(permuted)
  )
}

bootstrap_ci <- function(values, condition, seed) {
  nc_values <- values[condition == "NC"]
  pd_values <- values[condition == "PD"]
  set.seed(seed)
  deltas <- replicate(
    bootstrap_repetitions,
    mean(sample(pd_values, length(pd_values), replace = TRUE)) -
      mean(sample(nc_values, length(nc_values), replace = TRUE))
  )
  unname(quantile(deltas, c(0.025, 0.975), names = FALSE, type = 6))
}

summarise_test <- function(data, value_column, seed_offset = 0L) {
  value <- data[[value_column]]
  condition <- data$condition
  test <- exact_permutation_test(value, condition)
  ci <- bootstrap_ci(value, condition, base_seed + seed_offset)
  data.frame(
    n_NC = sum(condition == "NC"),
    n_PD = sum(condition == "PD"),
    mean_NC = mean(value[condition == "NC"]),
    mean_PD = mean(value[condition == "PD"]),
    median_NC = median(value[condition == "NC"]),
    median_PD = median(value[condition == "PD"]),
    delta_mean_PD_minus_NC = test$delta_mean,
    bootstrap_95CI_low = ci[[1]],
    bootstrap_95CI_high = ci[[2]],
    exact_p_two_sided = test$p_two_sided,
    exact_p_one_sided_PD_higher = test$p_one_sided_PD_higher,
    n_exact_permutations = test$n_permutations
  )
}

pathway_statistics <- bind_rows(lapply(seq_along(pathways_tested), function(i) {
  pathway_i <- pathways_tested[[i]]
  data_i <- pathway_scores %>% filter(as.character(pathway_name) == pathway_i)
  cbind(
    data.frame(pathway_name = pathway_i),
    summarise_test(data_i, "pathway_score", seed_offset = i)
  )
})) %>%
  mutate(
    BH_q_two_sided = p.adjust(exact_p_two_sided, method = "BH"),
    BH_q_one_sided_PD_higher = p.adjust(exact_p_one_sided_PD_higher, method = "BH")
  ) %>%
  arrange(match(pathway_name, pathways_tested))

write.csv(
  pathway_statistics,
  file.path(out_dir, "05_pathway_exact_permutation_statistics.csv"),
  row.names = FALSE
)

interaction_names <- unique(sample_scores$interaction_name)
lr_statistics <- bind_rows(lapply(seq_along(interaction_names), function(i) {
  interaction_i <- interaction_names[[i]]
  data_i <- sample_scores %>% filter(interaction_name == interaction_i)
  identity <- data_i[1, c("interaction_name", "pathway_name", "ligand", "receptor")]
  cbind(identity, summarise_test(data_i, "probability_core", seed_offset = 100L + i))
})) %>%
  mutate(
    BH_q_two_sided_all_LR = p.adjust(exact_p_two_sided, method = "BH"),
    BH_q_one_sided_all_LR = p.adjust(exact_p_one_sided_PD_higher, method = "BH")
  ) %>%
  arrange(BH_q_two_sided_all_LR, desc(delta_mean_PD_minus_NC))

write.csv(
  lr_statistics,
  file.path(out_dir, "06_LR_exact_permutation_statistics.csv"),
  row.names = FALSE
)

format_probability <- function(value) {
  ifelse(value < 0.001, "< 0.001", sprintf("= %.3f", value))
}

plot_labels <- pathway_statistics %>%
  transmute(
    pathway_name = factor(pathway_name, levels = pathways_tested),
    label = paste0(
      "exact P ", format_probability(exact_p_two_sided),
      "\nBH q ", format_probability(BH_q_two_sided)
    )
  ) %>%
  left_join(
    pathway_scores %>%
      group_by(pathway_name) %>%
      summarise(y = max(pathway_score) * 1.14 + 1e-6, .groups = "drop"),
    by = "pathway_name"
  )

group_colors <- c(NC = "#00BFC4", PD = "#F8766D")
p <- ggplot(pathway_scores, aes(x = condition, y = pathway_score, fill = condition)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.72) +
  geom_jitter(aes(color = condition), width = 0.10, size = 2.2, alpha = 0.90) +
  geom_text(
    data = plot_labels,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    lineheight = 0.95
  ) +
  facet_wrap(~pathway_name, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "ECM_Fib to Activated_Fib: sample-level differential statistics",
    subtitle = sprintf(
      "Independent samples: NC = %d, PD = %d; all eligible samples from %s",
      sum(eligible$condition == "NC"),
      sum(eligible$condition == "PD"),
      unique(eligible$dataset)
    ),
    x = NULL,
    y = "CellChat-compatible core pathway score"
  )

ggsave(
  file.path(out_dir, "Fig8F_ECM_to_Activated_sample_level_statistics.pdf"),
  p, width = 10.5, height = 4.2, useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "Fig8F_ECM_to_Activated_sample_level_statistics.png"),
  p, width = 10.5, height = 4.2, dpi = 600
)

significant_pd_table <- pathway_statistics %>%
  filter(BH_q_two_sided < 0.05, delta_mean_PD_minus_NC > 0) %>%
  select(pathway_name, BH_q_two_sided)

result_sentence <- if (nrow(significant_pd_table) > 0L) {
  significant_text <- paste0(
    significant_pd_table$pathway_name,
    " (BH q ", format_probability(significant_pd_table$BH_q_two_sided), ")"
  )
  non_significant_table <- pathway_statistics %>%
    filter(!(pathway_name %in% significant_pd_table$pathway_name))
  non_significant_text <- if (nrow(non_significant_table) > 0L) {
    paste0(
      " ", paste(non_significant_table$pathway_name, collapse = ", "),
      " showed the same direction but did not meet the two-sided BH-adjusted threshold",
      " (BH q ",
      paste(format_probability(non_significant_table$BH_q_two_sided), collapse = ", "),
      ")."
    )
  } else {
    ""
  }
  paste0(
    "At the biological-sample level, the ",
    paste(significant_text, collapse = ", "),
    " pathway scores were higher in PD after two-sided exact permutation testing and BH correction.",
    non_significant_text
  )
} else {
  paste0(
    "None of the four prespecified ECM pathway scores remained significant after ",
    "two-sided exact permutation testing and BH correction at the biological-sample level."
  )
}

report <- c(
  "# Figure 8F sample-level differential statistics",
  "",
  sprintf(
    "Eligible biological samples: NC = %d and PD = %d; each contained at least %d %s cells and %d %s cells. All eligible samples were from %s.",
    sum(eligible$condition == "NC"), sum(eligible$condition == "PD"),
    min_cells, sender, min_cells, receiver, unique(eligible$dataset)
  ),
  "",
  result_sentence,
  "",
  "## Figure 8F legend",
  "",
  paste0(
    "Figure 8F. Sample-level differential analysis of ECM_Fib-to-Activated_Fib ligand-receptor availability in NC and PD. ",
    "Each point represents an independent biological sample containing at least 10 cells of both populations (NC, n = ",
    sum(eligible$condition == "NC"), "; PD, n = ", sum(eligible$condition == "PD"),
    "; all from ", unique(eligible$dataset), "). For each displayed CellChat interaction, ligand expression in ECM_Fib and receptor or receptor-complex expression in Activated_Fib were summarized within each sample using the trimean and converted to the core CellChat mass-action probability using the default Hill parameters (Kh = 0.5, n = 1), without population-size weighting. Pathway scores are sums of the corresponding interaction probabilities. Boxes show medians and interquartile ranges. Two-sided P values were obtained by exhaustive permutation of the NC/PD labels across all 5,005 possible assignments; q values were adjusted across the four prespecified pathways using the Benjamini-Hochberg method. This analysis tests sample-level expression-based ligand-receptor availability and does not demonstrate physical signaling or causal fibroblast activation."
  ),
  "",
  "## Methods text",
  "",
  paste0(
    "To provide biological-replicate inference for the ECM_Fib-to-Activated_Fib axis, we performed a sample-level sensitivity analysis using the same composition-matched CellChat input used for Figure 8. Samples with at least 10 ECM_Fib and 10 Activated_Fib cells were retained (NC, n = ",
    sum(eligible$condition == "NC"), "; PD, n = ", sum(eligible$condition == "PD"),
    "; all from ", unique(eligible$dataset), "). Ligand and receptor expression was summarized separately for each cell population and sample using the trimean. Multi-subunit receptor expression was represented by the geometric mean of its detected subunits. The ligand-receptor product was transformed using the core CellChat Hill function (Kh = 0.5, n = 1), without population-size weighting. Probabilities were summed within the prespecified COLLAGEN, FN1, LAMININ, and PERIOSTIN pathways. NC-PD differences were assessed using an exhaustive two-sided sample-label permutation test (5,005 unique assignments), with Benjamini-Hochberg correction across the four pathways. Stratified nonparametric bootstrap resampling of biological samples (10,000 repetitions) was used to obtain 95% confidence intervals for the difference in group means."
  ),
  "",
  "## Discussion text",
  "",
  paste0(
    "The CellChat-identified enhancement of ECM_Fib-to-Activated_Fib communication in PD (Figure 8E-F) is an expression-based inference derived from matched cell numbers without population-size weighting. ",
    result_sentence, " The inferred communication probabilities reflect ligand-receptor transcript availability under the CellChat mass-action model rather than direct evidence of ligand secretion, receptor activation, physical cell proximity, or causal conversion of ECM_Fib into Activated_Fib. These structured predictions prioritize experimental tests of ECM-mediated fibroblast activation, including co-culture perturbation of CD44 and ITGAV/ITGB5 signaling."
  ),
  "",
  "## Interpretation boundary",
  "",
  "The aggregate CellChat objects and the sample-level analysis answer different questions. The aggregate objects estimate pathway structure after cell-number matching; the exact permutation test evaluates whether the corresponding expression-based availability score differs across independent samples. Neither analysis establishes molecular binding or causality."
)
writeLines(report, file.path(out_dir, "07_manuscript_text_and_interpretation.md"))

writeLines(
  capture.output({
    cat("CellChat version:", as.character(packageVersion("CellChat")), "\n")
    cat("dplyr version:", as.character(packageVersion("dplyr")), "\n")
    cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n")
    cat("Minimum cells per population per sample:", min_cells, "\n")
    cat("Eligible NC samples:", sum(eligible$condition == "NC"), "\n")
    cat("Eligible PD samples:", sum(eligible$condition == "PD"), "\n")
    cat("Exact label assignments:", choose(nrow(eligible), sum(eligible$condition == "PD")), "\n")
    cat("Bootstrap repetitions:", bootstrap_repetitions, "\n")
    cat("Base seed:", base_seed, "\n")
    print(sessionInfo())
  }),
  file.path(out_dir, "08_sessionInfo.txt")
)

message("Completed: ", normalizePath(out_dir))
print(pathway_statistics)
