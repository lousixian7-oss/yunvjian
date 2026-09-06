# Prediction reconciliation — 2026-09-05

Repository: https://github.com/lousixian7-oss/yunvjian

Figure 5E and Figure S2E/G now use the Table S6 preprocessing path and the same fixed probabilities. The saved logistic model was not refitted or recalibrated. All four cohorts reproduce the S6 AUC, Brier score, calibration intercept and slope within 1e-8. Confidence intervals are retained from the archived bootstrap results, not recomputed. GSE223924 predictions are read from the archived MOR results; this runner does not repeat RNA-seq normalization.

The GSE16134 cohort uses its own gene means and sample SDs with the algebraically equivalent training-fitted standardized coefficients. This is cohort-based transportability evaluation; it does not establish a frozen single-new-patient deployment pipeline. GSE16134 AUC is 0.9225449516 and Brier score is 0.1443679409. The historical raw-expression figure path had Brier 0.0974642744 and is superseded. The S6 path was selected for consistency with the declared model specification, not because it improves performance.

## Reproduce

From the repository root, with R and ggplot2 installed:

    Rscript figure5/probability_reconciliation/reconcile_probabilities.R . figure5/probability_reconciliation

Then, with Python packages pillow, pypdf and reportlab, and Poppler pdftoppm on PATH (or PDFTOPPM pointing to its executable):

    python figure5/probability_reconciliation/assemble_reconciled.py

The shipped panels preserve the latest composite layout and unchanged panels. The R runner rebuilds only the three affected panels. Use this runner and assembler for the reconciled release; the historical whole-figure script retains its older assembly workflow. Local package versions are recorded in R_sessionInfo.txt. A clean-machine end-to-end rerun of every manuscript analysis has not been performed.

Unified_fixed_predictions.csv is the common prediction source; Verified_S6_metrics.csv is the numerical check. Table S6 numerical cells are unchanged. Upload all files from the patch preserving repository-relative paths; do not place github_patch itself inside the repository. This release includes the reconciled fixed predictions and corrected overlapping-cohort labels.

## Overlap interpretation

Full GSE16134 is an overlapping-cohort comparison. The 66-site/30-participant post hoc subset and Figure S10 are documented in [the subset analysis](../nonoverlap_validation_20260905/README.md). These are nested analyses, not two independent external cohorts.
