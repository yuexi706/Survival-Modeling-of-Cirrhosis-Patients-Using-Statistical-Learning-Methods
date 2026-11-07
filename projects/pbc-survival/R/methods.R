# Endpoint: time to death, censoring at transplantation or administrative follow-up.
# This is cause-specific discrimination, not calibrated cumulative mortality risk.
pbc_numeric <- c("age", "bili", "chol", "albumin", "copper", "alk.phos",
                 "ast", "trig", "platelet", "protime")
pbc_levels <- list(sex = c("f", "m"), ascites = c("0", "1"),
                   hepato = c("0", "1"), spiders = c("0", "1"),
                   edema = c("0", "0.5", "1"), stage = as.character(1:4))
pbc_predictors <- c(pbc_numeric, names(pbc_levels))
pbc_logs <- c("bili", "chol", "copper", "alk.phos", "ast", "trig", "protime")
cox_formula <- survival::Surv(time, death) ~ age + sex + bili + albumin + protime + edema

fit_pbc_preprocessor <- function(train) {
  medians <- vapply(train[pbc_numeric], stats::median, numeric(1), na.rm = TRUE)
  modes <- vapply(names(pbc_levels), function(v) {
    z <- table(factor(train[[v]], levels = pbc_levels[[v]]))
    if (!sum(z)) stop("No observed training values: ", v)
    names(z)[which.max(z)]
  }, character(1))
  stopifnot(all(is.finite(medians)))
  list(medians = medians, modes = modes)
}

apply_pbc_preprocessor <- function(d, prep) {
  for (v in pbc_numeric) {
    d[[v]][is.na(d[[v]])] <- prep$medians[[v]]
    if (v %in% pbc_logs) d[[v]] <- log1p(d[[v]])
  }
  for (v in names(pbc_levels)) {
    z <- as.character(d[[v]])
    if (any(!is.na(z) & !z %in% pbc_levels[[v]])) stop("Unexpected category: ", v)
    z[is.na(z)] <- prep$modes[[v]]
    d[[v]] <- factor(z, levels = pbc_levels[[v]])
  }
  d$death <- as.integer(d$status == 2L)
  d
}

pbc_matrix <- function(d) {
  # Explicit allowlist: id, time, status and treatment assignment are never predictors.
  stats::model.matrix(~ ., data = d[pbc_predictors])[, -1, drop = FALSE]
}

risk_cindex <- function(time, death, risk) {
  if (!sum(death) || length(unique(risk)) < 2L) return(0.5)
  # Greater risk means shorter survival. The explicit reversal is essential.
  as.numeric(survival::concordance(survival::Surv(time, death) ~ risk,
                                  reverse = TRUE)$concordance)
}

tune_lasso <- function(raw_train, seed, k = 4L) {
  folds <- stratified_folds(raw_train$status, k, seed)
  lambda <- exp(seq(log(1), log(0.001), length.out = 40))
  scores <- matrix(NA_real_, k, length(lambda))
  for (j in seq_len(k)) {
    raw_a <- raw_train[folds != j, ]; raw_b <- raw_train[folds == j, ]
    prep <- fit_pbc_preprocessor(raw_a)
    a <- apply_pbc_preprocessor(raw_a, prep); b <- apply_pbc_preprocessor(raw_b, prep)
    fit <- glmnet::glmnet(pbc_matrix(a), survival::Surv(a$time, a$death),
                          family = "cox", alpha = 1, lambda = lambda)
    risk <- as.matrix(predict(fit, newx = pbc_matrix(b), s = lambda, type = "link"))
    scores[j, ] <- apply(risk, 2, function(r) risk_cindex(b$time, b$death, r))
  }
  # Preprocessing is re-estimated inside every inner fold, including imputation.
  means <- colMeans(scores)
  best <- which.max(means)
  se <- stats::sd(scores[, best]) / sqrt(k)
  one_se <- which(means >= means[best] - se)[1L] # lambda sorted strongest first
  list(lambda_min = lambda[best], lambda_1se = lambda[one_se],
       curve = data.frame(lambda, mean_cindex = means,
                          se = apply(scores, 2, stats::sd) / sqrt(k)))
}
