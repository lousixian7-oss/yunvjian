# Yunvjian periodontitis analysis archive

Analysis scripts and processed core data supporting the current Yunvjian manuscript (2026-09-05). Main figures use `figure1` to `figure9`; supplementary figures use `figures1` to `figures8`. Each folder explains its source workflow and limitations.

## Data sources

Bulk discovery: GSE23586, GSE156993, GSE10334. Additional evaluation: full GSE16134 (overlapping-cohort comparison), its post hoc participant-disjoint subset, and exploratory GSE223924. Single-cell RNA-seq: GSE152042, GSE164241, GSE171213. Spatial transcriptomics: GSE206621 (GSM6258255–GSM6258258).

The full processed discovery matrix has 172 samples; WGCNA excludes GSM261278 and uses 171 samples. The WGCNA network then uses the 5,000 most variable eligible genes. Do not use the 171-sample number for the full discovery dataset.

## How to use

Start with the README in the relevant figure folder. `provenance.csv` maps archived files to original local sources with SHA-256 checksums. `input_references.csv` lists input and working-directory statements recovered from scripts; these are references, not proof that all dependencies are included. `R_packages_referenced.txt` lists detected package dependencies. Existing scripts retain original code, relative-path assumptions, and historical paths; configure those paths and use separate output directories before rerunning.

Core data include compressed expression matrices, sample annotations, cell metadata, fibroblast/hdWGCNA/model/CellChat result tables and spatial RCTD weights. Large Seurat objects are represented by metadata and embedding exports plus source records, as permitted by the requested data checklist. Metadata cannot replace raw counts for differential expression or rerunning integration.

## Reproducibility status

This release organizes existing scripts and saved results; it does not claim that all analyses have been rerun in a clean environment. The final three-gene model is AGT/CXCR4/FOS. Historical alternative models and exploratory screening must not be interpreted as independent confirmatory selection. Figure 1 was created in BioRender. Cytoscape/STRING/BATMAN steps were performed through graphical interfaces; no scripts are required for those manual steps. The final global kBET scripts, logs and result tables have been recovered from the original analysis folder in the Recycle Bin and are archived under figures3 and coredata/batch_integration; see figures3/README.md for run order and provenance.

See `RELEASE_CHECKS.json` for file counts, data dimensions, checksums and remaining gaps. No software license is asserted for third-party code or upstream datasets; their original terms remain applicable.


## 2026-09-05 prediction reconciliation

Figure 5E and Figure S2E/G now share the fixed Table S6 predictions. See [reconciliation instructions](figure5/probability_reconciliation/README.md). Use the dated reconciled figures and runner. No model refitting or recalibration was performed.

## 2026-09-06 overlap audit and reporting update

The full GSE16134 series contains 244 sites from 90 patients represented in development. Its 310-site result is retained for transparency as an overlapping-cohort comparison. Excluding shared patients leaves 66 sites from 30 patients (60 affected, six unaffected); this same-source post hoc subset is not an independently recruited or previously unseen validation cohort. It is not added again to the sample total.

[Figure S10 and the complete analysis snapshot](figure5/nonoverlap_validation_20260905/README.md) provide predictions, preprocessing comparisons, patient-cluster bootstrap results, raw-file hashes, and overlap disposition. Primary AUC is 0.986 (95% CI 0.950-1.000), but Brier score is 0.205 and calibration is poor. Only six unaffected sites are available; calibration-slope bootstrap fits are unstable. Model coefficients were not refitted or recalibrated. Figure 5E now labels the overlap. The original probability-reconciliation update is included in this release.

Updated Tables S1/S6 are in [reporting tables](figure5/nonoverlap_validation_20260905/reporting_tables). The local manuscript and TRIPOD were synchronized separately; their unpublished DOCX files are not included in this code archive.
