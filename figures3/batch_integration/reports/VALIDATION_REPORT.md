## Material Passport

- Origin Skill: experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-21
- Verification Status: VERIFIED
- Version Label: batch_integration_validation_v1

## Validation Report

- Source: local Seurat object and reproducible R scripts in this directory
- Overall Confidence: CAUTION

### Statistical findings

| Metric | Value | Interpretation | Confidence |
|---|---:|---|---|
| Official kBET, PCA before | 0.940 | Strong local dataset structure | SOLID |
| Official kBET, explicit Harmony | 0.571 | Substantial improvement, incomplete removal | CAUTION |
| iLISI, PCA to explicit Harmony | 1.532 to 2.212 | Increased local dataset diversity | SOLID |
| Cross-dataset 30-NN | 0.268 to 0.509 | Increased cross-dataset mixing | SOLID |
| ARI, original vs explicit clusters | 0.952 | Strong clustering stability | SOLID |
| Sample cluster-profile correlation | median 0.999 | Sample composition is highly stable | SOLID |

### Warnings

| Type | Detail | Affected results |
|---|---|---|
| Residual batch structure | Official kBET remains above its expected rejection rate | Especially pericytes and macrophages |
| Composition limitation | At least one expected local batch count is below 5 | Cycling cells, mast cells, pDCs |
| Sensitivity of kBET | Large cell numbers make small local deviations detectable | Interpretation of absolute rejection rate |
| Pseudoreplication | Cells are not independent biological replicates | NC versus PD differential analysis |
| QC reconstruction | Only post-filter object is available | Per-sample QC pass/removal rates cannot be recovered |

### Fallacy scan

- Coverage: 11/11 fallacy types checked.
- Simpson's paradox: CAUTION. Overall batch summaries can obscure cell-type-specific behavior; addressed by stratification.
- Ecological fallacy: CAUTION. Dataset-level mixing cannot establish cell-level disease effects.
- Berkson's paradox: NOTE. The object contains only QC-passing cells; raw exclusion effects cannot be reconstructed.
- Collider bias: NOTE. No causal adjustment claim was made; QC conditioning should still be described.
- Base-rate neglect: NOTE. Not a diagnostic-accuracy analysis; cell-type dataset proportions were reported.
- Regression to the mean: NOTE. Not a pre-post intervention analysis.
- Survivorship bias: CAUTION. Post-QC-only availability prevents direct attrition assessment.
- Look-elsewhere effect: CAUTION. Eleven cell types were examined; results are descriptive and not selective significance claims.
- Garden of forking paths: CAUTION. Integration parameters and random seeds are now explicitly recorded; sensitivity analysis is provided.
- Correlation is not causation: CAUTION. Batch-mixing and disease associations must not be phrased causally.
- Reverse causality: NOTE. Not directly applicable to integration diagnostics; disease associations remain observational.

### Reproducibility

- Method: fixed-seed rerun with explicit Harmony, graph reconstruction, formal kBET 0.99.6, and file-level validation.
- Verdict: REPRODUCIBLE.
- The official kBET result closely agreed with the independently implemented kBET-style sensitivity result: PCA 0.940 versus 0.934 and explicit Harmony 0.571 versus 0.561.

