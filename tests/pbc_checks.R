require_packages("survival")
source("projects/pbc-survival/R/methods.R")
raw <- survival::pbc
prep <- fit_pbc_preprocessor(raw[1:300, ])
holdout <- raw[301:418, ]; holdout$bili <- 1e12
stopifnot(identical(prep, fit_pbc_preprocessor(raw[1:300, ])))
# Missing held-out values must receive the TRAINING median (before log1p).
holdout$bili[1] <- NA
processed <- apply_pbc_preprocessor(holdout, prep)
stopifnot(processed$bili[1] == log1p(median(raw$bili[1:300], na.rm = TRUE)))
x <- pbc_matrix(processed)
stopifnot(!any(grepl("^(id|time|status|death|trt)$", colnames(x))))
stopifnot(identical(processed$death, as.integer(holdout$status == 2)))
stopifnot(risk_cindex(1:6, rep(1, 6), 6:1) == 1,
          risk_cindex(1:6, rep(1, 6), 1:6) == 0)

if ("--results" %in% commandArgs(trailingOnly=TRUE)) {
 pbc_dir <- "projects/pbc-survival/results"
 pred <- read.csv(file.path(pbc_dir, "out_of_fold_predictions.csv"))
 # One out-of-fold score per patient/model/repeat; no training predictions mixed in.
 stopifnot(!anyDuplicated(pred[c("repeat_id", "id", "model")]))
 counts <- aggregate(id ~ repeat_id + model, pred, length)
 stopifnot(all(counts$id == 418))
 metrics <- read.csv(file.path(pbc_dir, "fold_metrics.csv"))
 stopifnot(all(is.finite(metrics$cindex)), all(metrics$cindex >= 0 & metrics$cindex <= 1))
 ci <- read.csv(file.path(pbc_dir, "five_year_estimands.csv"))
 stopifnot(ci$death_cif + ci$transplant_cif <= 1,
           ci$death_cif <= ci$naive_one_minus_km)

}
