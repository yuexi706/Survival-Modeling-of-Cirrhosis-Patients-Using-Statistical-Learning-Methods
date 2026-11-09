# Methodological audit

Based on supplied rendered code and reports. The original CSV is missing; original results are not relabeled as extension results.

## PBC survival: selected as the main extension

The original code fits median/mode imputation to the full data before splitting. Its LASSO feature construction excludes follow-up and outcomes but does not exclude the `ID` field. Its forest and Cox formulas explicitly name clinical predictors, so the ID issue applies specifically to the LASSO design. The code creates risk tertiles across all records after training on a subset, mixing development and evaluation samples in risk-group displays. A single split is insufficient to establish a small performance advantage. Transplantation is treated as censoring for cause-specific evaluation and as a competing event in a separate Fine-Gray section; these targets need distinct interpretation.

The extension fits preprocessing within each resampling level, removes ID by allowlist, uses repeated nested evaluation for LASSO and measures selection stability. It shows cumulative incidence separately from cause-specific discrimination. Fine-Gray predictive validation was not added, and no calibrated competing-risk prediction claim is made.

