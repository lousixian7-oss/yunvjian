# Figure S3 — Global single-cell QC and integration

The final original kBET workflow has been recovered and archived in `batch_integration/`. Its original directory was absent from the active project because it was in the G-drive Recycle Bin. The recovered `03_official_kbet.R` exactly matches the successfully applied historical creation and final correction (`n_repeat <- 100`). This release recovers existing code and results; it does not rerun the analyses.

## Run order

1. `scripts/01_batch_integration_revision.R`: original PCA -> explicit Harmony and comparison embeddings.
2. `scripts/02_cluster_stability_and_kbet.R`: clustering stability and the separate custom kBET-style diagnostic.
3. `scripts/03_official_kbet.R`: official kBET v0.99.6 evaluation and the final Figure S3F panel.

Scripts take input RDS, output directory and an existing R-library path as command-line arguments. The initial input was `G:/1Yunvjian/0a26.7.7singlecell/8.10/yjsl_clean_allGSE.rds`, which still exists locally. Script 01 produces the object containing `harmony_explicit`, required by script 03. The earlier global metadata/embedding exports in coredata represent the original integration and do not substitute for that intermediate object.

Example after configuring paths:

```text
Rscript scripts/01_batch_integration_revision.R <original_rds> <output_dir> <existing_R_library>
Rscript scripts/02_cluster_stability_and_kbet.R <integrated_rds> <output_dir> <existing_R_library>
Rscript scripts/03_official_kbet.R <integrated_rds> <output_dir> <existing_R_library>
```

`00_install_packages.R` and `historical_launchers/` are retained as historical source records, not the recommended way to install or launch this archive. Reuse an existing compatible library; adjust paths rather than blindly running the old launchers. No packages were installed during recovery.

## Matched final results

Official weighted rejection rates: PCA 0.9404629877; existing Harmony 0.6028385616; explicit Harmony 0.5708235463. Lower rejection indicates better local dataset mixing. This tests mixing within cell types, rather than differences between disease groups. It supports substantial improvement with residual dataset structure.

The official run used PCA dimensions 1–15, a maximum of 500 cells per cell-type/dataset stratum, 30 neighbors, seed 20260821, and 100 repetitions. The log records successful completion on 2026-08-21 at 12:48:47 CST. Saved official and custom diagnostic results are distinguished in `coredata/batch_integration/`. Environment records and the original final PDF are included here.
