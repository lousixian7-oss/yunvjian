# Batch integration revision: results summary

## Verification status

- Status: completed and re-opened successfully
- Input cells: 118,050
- Input genes: 54,492
- Software: R 4.5.2, Seurat 5.5.1, SeuratObject 5.4.0, Harmony 2.0.5
- Harmony covariate: `dataset`
- Dimensions: PCA 1-15
- Explicit Harmony maximum iterations: 10
- Observed convergence: 5 iterations
- Random seed: 42

## Overall batch-mixing diagnostics

Metrics were calculated on a fixed-seed, cell-type-stratified sample capped at 500 cells per dataset-cell-type stratum. These diagnostics are not presented as formal kBET results.

| Reduction | Cross-dataset fraction among 30-NN | Cross/expected ratio | iLISI | Dataset silhouette |
|---|---:|---:|---:|---:|
| PCA before integration | 0.2679 | 0.4243 | 1.5324 | 0.0462 |
| Existing Harmony | 0.5033 | 0.7920 | 2.1649 | -0.0271 |
| Explicit Harmony rerun | 0.5094 | 0.7975 | 2.2117 | -0.0339 |

The explicit Harmony rerun improved all three mixing summaries relative to PCA and was slightly better than the existing Harmony embedding. The existing and explicit Harmony dimensions were highly similar (mean dimension-wise Pearson correlation 0.9832; range 0.9661-0.9991).

## Cell-type-specific interpretation

- Best mixing among major lineages: neutrophils (cross/expected ratio 0.910), mast cells (0.863), epithelial cells (0.837), fibroblasts (0.826), and endothelial cells (0.825).
- Pericytes retained the strongest local dataset structure (cross/expected ratio 0.607; dataset silhouette 0.0056).
- Cycling cells had low absolute mixing because GSE171213 contained only 22 sampled cycling cells; this metric is composition-limited and should not be interpreted as batch failure by itself.
- The prominent GSE164241 regions in the dataset-colored UMAP largely coincide with cell-type composition imbalance: GSE164241 contributed 95.3% of cycling cells, 87.8% of plasma cells, 87.3% of mast cells, 86.9% of fibroblasts, and 84.1% of B cells.

## Decision

The original reviewer observation that dataset enrichment is visible is correct, but the new analysis does not support the stronger conclusion that Harmony failed globally. Dataset mixing improved substantially after Harmony, while some lineage-specific residual structure remains, especially in pericytes. The recommended manuscript strategy is to retain the explicit dataset-only Harmony result as the main integration, report the quantitative metrics, and use cell-type-stratified or sample-level sensitivity analyses rather than switching directly to scVI.

## Official kBET validation

The official kBET 0.99.6 implementation was run within cell-type strata on the same fixed-seed balanced sample, using 30 nearest neighbors and 100 repeats.

| Reduction | Weighted observed rejection rate | Weighted expected rejection rate |
|---|---:|---:|
| PCA before integration | 0.9405 | 0.0454 |
| Existing Harmony | 0.6028 | 0.0429 |
| Explicit Harmony rerun | 0.5708 | 0.0443 |

Harmony substantially reduced kBET rejection relative to PCA, while the remaining rejection rate indicates incomplete batch removal. The official and independently implemented kBET-style results agreed closely.

## Clustering stability after explicit Harmony graph reconstruction

- Original clusters: 17; explicit Harmony clusters: 17.
- Adjusted Rand index: 0.9517.
- Normalized mutual information: 0.9435.
- Cell mapping accuracy after majority mapping: 0.9629.
- New-cluster cell-type purity: 0.9634.
- Median sample-level cluster-profile correlation: 0.9986.
- Median sample-level total variation distance: 0.0208.

All four prespecified stability rules passed. The explicit Harmony rerun therefore supports the robustness of the original biological structure and does not, by itself, require rebuilding every downstream analysis.

## Important limitation

Only post-filter cells were available in the saved object. Raw per-sample QC pass rates and the percentage removed for `percent.mt > 30` cannot be reconstructed from this object and require the pre-filter object or raw matrices.
