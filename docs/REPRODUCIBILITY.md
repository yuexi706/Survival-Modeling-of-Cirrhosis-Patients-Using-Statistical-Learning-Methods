# Reproducibility record

The extension was executed locally on R 4.4.1, macOS arm64. The committed sessionInfo.txt and renv.lock record package versions. A clean-machine restore has not been verified.

```sh
Rscript scripts/install.R
Rscript tests/validate.R
Rscript projects/pbc-survival/run.R
Rscript tests/validate.R --results
python3 scripts/build_reports.py
```

Three repetitions, five outer folds, four LASSO inner folds. Each relevant training fold gets its own imputation fit. Lambda grid: 40 log-spaced values from 1 to 0.001; forest: 500 trees, mtry=4, minimum node size=15, one thread. Seeds, feature schema and transformations are explicit. PBC_REPEATS changes the experiment; committed results use the default 3.

Tests check held-out imputation, predictor exclusion, event mapping, concordance direction, one out-of-fold prediction per record/model/repetition and valid competing-event probabilities. Fold C-indices are computed separately; overlapping folds produce descriptive variation, not independent inferential intervals.

The GitHub workflow runs fast statistical invariants only. It does not refit models, certify conclusions or test exact lockfile restoration. Full model execution and result validation were done locally. Small numerical differences can occur across systems despite seeds.
