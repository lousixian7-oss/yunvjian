# Corrected three-gene figures, 2026-09-02

Use Figure5_three_gene_corrected_20260902 and FigureS2_three_gene_corrected_20260902 (PDF, PNG, TIFF). Individual vector PDF panels are in panels/. Old releases are preserved. Word documents were not changed.

## What changed

- S2 A-C use the same saved sample predictions from a single 41-gene, eight-algorithm run, not mixed 40/41-gene sources. NR3C1 is the additional candidate.
- The seed-123 sample-stratified 70/30 split exactly matches the original 8.12 sample indices: 122 training and 50 screening-holdout samples. All eight training and holdout AUCs reproduce the recorded 41-gene workbook values.
- This screening split differs from the patient-grouped split used by the downstream three-gene logistic model. The screening holdout is also used for permutation importance, so it is not independent validation of feature selection.
- The model previously called LASSO selected alpha=0.55 and lambda=0.0658314 using caret's glmnet tuning. The figure now correctly labels it Elastic net; no forced alpha=1 refit was substituted.
- S2 D shows the actual rerun ranks for AGT, CXCR4, FOS and comparator EPHX2. It replaces the previously asserted stable four-method intersection.
- The rerun four-method Top-20 overlap is CXCR4, FOS, SLC6A4; the earlier workbook reported AGT, CXCR4, FOS, EPHX2. AUC reproducibility does not establish permutation-ranking reproducibility. Permutation random state and tie handling affect rankings, and the historical permutation random state was not saved. No seeds were searched to obtain a preferred signature.
- Ranks use decreasing mean DALEX permutation dropout RMSE. Ties are resolved by candidate-list order, not interpreted as different importance.
- S2 F was recovered from the original WGCNA .RData (171 samples and seven eigengenes). Its CXCR4/FOS values and module labels match the original heatmap after rounding. AGT was absent from the 5,000-gene WGCNA input; absence does not mean zero correlation.
- S2 E/G retain the previously calculated three-gene probabilities and DCA. Figure 5 numerical results are retained. Its nomogram labels/tick density were repaired with the saved model and verified against fitted probabilities.

The user-selected AGT/CXCR4/FOS model is retained as a reduced candidate model, not relabelled as a newly established or unchanged algorithmic intersection. Retrospective evaluation of several signatures is not prospective validation. Figure S2 D and the Elastic net terminology require matching captions wherever the new figure is used.

## Saved evidence

- candidate41.txt: exact CytoNCA candidate order.
- rerun_41gene.R / rerun.log: reproducible original-method fitting pipeline.
- models/: caret fits and DALEX permutation objects.
- eight_models_and_data.rds: all fits, input tables, genes and post-fit RNG state.
- eight_models_sample_predictions.csv: all 172 sample predictions for each model.
- eight_models_metrics.csv / best_tuning_parameters.txt: numerical metrics and selected tuning.
- eight_models_importance_rankings.csv / *_importance_raw.csv: full ranking audit.
- WGCNA_correlations_recovered.csv: source-derived Pearson coefficients.
- plot_corrected_S2.R / fix_nomogram.R / assemble_figures.py: figure builders.
- R_package_versions.csv / R_sessionInfo.txt: actual environment.
- QA.md / figure_export_checks.json: final export verification.

The required packages were already usable in the existing user library when accessed outside the restricted runtime. No package download, reinstall, or update was performed.
