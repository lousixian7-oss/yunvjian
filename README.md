# Yunvjian periodontitis analysis archive

Analysis scripts and processed core data supporting the current Yunvjian manuscript (2026-09-05). Main figures use `figure1` to `figure9`; supplementary figures use `figures1` to `figures8`. Each folder explains its source workflow and limitations.

## Data sources

Bulk discovery: GSE23586, GSE156993, GSE10334. External evaluation: GSE16134 and GSE223924. Single-cell RNA-seq: GSE152042, GSE164241, GSE171213. Spatial transcriptomics: GSE206621 (GSM6258255–GSM6258258).

The full processed discovery matrix has 172 samples; WGCNA excludes GSM261278 and uses 171 samples. The WGCNA network then uses the 5,000 most variable eligible genes. Do not use the 171-sample number for the full discovery dataset.

## How to use

Start with the README in the relevant figure folder. `provenance.csv` maps archived files to original local sources with SHA-256 checksums. `input_references.csv` lists input and working-directory statements recovered from scripts; these are references, not proof that all dependencies are included. `R_packages_referenced.txt` lists detected package dependencies. Existing scripts retain original code, relative-path assumptions, and historical paths; configure those paths and use separate output directories before rerunning.

Core data include compressed expression matrices, sample annotations, cell metadata, fibroblast/hdWGCNA/model/CellChat result tables and spatial RCTD weights. Large Seurat objects are represented by metadata and embedding exports plus source records, as permitted by the requested data checklist. Metadata cannot replace raw counts for differential expression or rerunning integration.

## Reproducibility status

This release organizes existing scripts and saved results; it does not claim that all analyses have been rerun in a clean environment. The final three-gene model is AGT/CXCR4/FOS. Historical alternative models and exploratory screening must not be interpreted as independent confirmatory selection. GUI-based network steps and the separate final global kBET script remain incompletely captured. No workflow code was invented for Figure 1.

See `RELEASE_CHECKS.json` for file counts, data dimensions, checksums and remaining gaps. No software license is asserted for third-party code or upstream datasets; their original terms remain applicable.
