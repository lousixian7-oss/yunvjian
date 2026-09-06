# figure5: Fixed AGT CXCR4 FOS model

compare_three_and_ephx2_models.R -> compare_models_gse223924.R -> rebuild_Figure5_FigureS2_three_gene_20260902.R; corrected_release/fix_nomogram.R and assemble_figures.py apply final display corrections. Final model has three genes; historical alternatives are sensitivity analyses.

Scripts preserve their original contents and may contain author-specific working directories. Inspect `setwd`, input paths and output paths before execution. The repository is an archive of existing analyses, not a newly rerun end-to-end pipeline. See the root provenance and input-reference manifests.


## 2026-09-05 prediction reconciliation

Figure 5E and Figure S2E/G now share the fixed Table S6 predictions. See [reconciliation instructions](probability_reconciliation/README.md). Use the dated reconciled figures and runner. No model refitting or recalibration was performed.

## 2026-09-05 participant-overlap audit and sensitivity evaluation

Full GSE16134 is an overlapping-cohort comparison, not independent external validation. Figure S10 separately evaluates the post hoc 66-site/30-patient participant-disjoint subset. See [analysis snapshot and scope](nonoverlap_validation_20260905/README.md). Original model coefficients and full-cohort numerical results are unchanged. The updated Figure 5E title identifies the overlap. This update includes the completed local analysis and its reporting limitations.
