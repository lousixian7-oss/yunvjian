# Participant-disjoint sensitivity evaluation, 2026-09-05

This is an archival snapshot of the completed post hoc analysis, not a new independently recruited validation cohort. Full GSE16134 (310 sites/120 participants) contains 244 sites from the 90 patients represented in development. Excluding these patients retains 66 sites/30 participants (60 affected, six unaffected sites). The six unaffected sites are from periodontitis patients.

The primary path separately RMA-processes the 66 CEL files, aggregates probes using the original mapping, applies subset gene means/sample SDs and the unchanged training-standardized AGT/CXCR4/FOS coefficients. No coefficient refitting, threshold optimization or recalibration was performed. RMA310/subset scaling and filtered original probabilities are retained as comparisons, not selected on performance.

Primary AUC is 0.986111 (patient-cluster bootstrap 95% CI 0.950-1.000), Brier 0.205029, calibration intercept 4.011359 and slope 1.726559. Only 1,289 of 2,000 bootstrap slope fits converged. The high AUC does not establish reliable clinical probabilities or clinical utility.

## Files and reproduction

- `nonoverlap_validation/`: sample predictions, expression/scaling parameters, all bootstrap results, influence analysis, ROC/calibration/DCA data, raw CEL hashes and R session information.
- `overlap_audit/`: accession crosswalk, disposition of all 310 samples, unresolved mappings, audit scripts and reports.
- `prepare_nonoverlap.py`, `nonoverlap_validation.R`, `report_nonoverlap.py`, `figure_S10.R`: original executed workflow scripts.

These scripts retain the original Windows workspace paths. Original working directory: `G:/1Yunvjian/机器学习/8.10/结果`. The relative analysis root is `tripod_20260905`; the existing model/core-data archive is `G:/1Yunvjian/1A玉女煎/code`. When reproducing elsewhere, configure these paths, restore the required archived development/full-cohort inputs, and obtain the 66 CEL files identified by the retained-sample list and CEL manifest. Raw CEL files and large probe-level caches are not included in this lightweight update. Use the package versions in sessionInfo.txt; R 4.5.2 was used. This archive does not claim a clean-environment end-to-end rerun of every manuscript analysis.

Figure S10 uses the fixed primary prediction CSV. Original full-cohort Figure 5/S2 statistics remain unchanged; the two accompanying plotting scripts only correct the overlapping-cohort labels. See Tables S1/S6 and the revised manuscript/TRIPOD for reporting scope. This update includes the completed local analysis and its reporting limitations.
