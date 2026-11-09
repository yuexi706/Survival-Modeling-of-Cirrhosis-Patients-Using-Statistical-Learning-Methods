source("R/common.R")
require_packages(c("survival", "glmnet", "ranger", "cmprsk", "ggplot2"))
source("projects/pbc-survival/R/methods.R")
out <- "projects/pbc-survival/results"
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
raw <- survival::pbc
stopifnot(nrow(raw) == 418L, !anyNA(raw[c("id", "time", "status")]),
          all(raw$status %in% 0:2), all(raw$time > 0), !anyDuplicated(raw$id))
save_table(data.frame(variable = names(raw), missing_n = colSums(is.na(raw)),
                      missing_fraction = colMeans(is.na(raw))), file.path(out, "missingness.csv"))
save_table(as.data.frame(table(status = factor(raw$status, levels = 0:2,
  labels = c("Censored", "Transplant", "Death")))), file.path(out, "endpoint_counts.csv"))

repeats <- as.integer(Sys.getenv("PBC_REPEATS", "3"))
outer_k <- 5L
predictions <- list(); fold_metrics <- list(); selections <- list(); tuning <- list()
for (r in seq_len(repeats)) {
  outer <- stratified_folds(raw$status, outer_k, 20261107L + r)
  for (f in seq_len(outer_k)) {
    message("PBC repeat ", r, "/", repeats, ", fold ", f, "/", outer_k)
    raw_train <- raw[outer != f, ]; raw_test <- raw[outer == f, ]
    tuning_fit <- tune_lasso(raw_train, 10000L * r + f)
    prep <- fit_pbc_preprocessor(raw_train)
    train <- apply_pbc_preprocessor(raw_train, prep)
    test <- apply_pbc_preprocessor(raw_test, prep)
    x_train <- pbc_matrix(train); x_test <- pbc_matrix(test)
    baseline <- survival::coxph(cox_formula, data = train, x = TRUE)
    lasso <- glmnet::glmnet(x_train, survival::Surv(train$time, train$death),
      family = "cox", alpha = 1,
      lambda = sort(unique(c(tuning_fit$lambda_min, tuning_fit$lambda_1se)), decreasing = TRUE))
    forest_data <- train[c("time", "death", pbc_predictors)]
    forest <- ranger::ranger(survival::Surv(time, death) ~ ., data = forest_data,
      num.trees = 500, mtry = 4, min.node.size = 15, seed = 90000L + 100L*r + f,
      num.threads = 1, importance = "permutation")
    forest_pred <- predict(forest, data = test[pbc_predictors], num.threads = 1)
    risks <- list(
      Cox = as.numeric(predict(baseline, newdata = test, type = "lp")),
      LASSO = as.numeric(predict(lasso, newx = x_test, s = tuning_fit$lambda_min)),
      LASSO_1SE = as.numeric(predict(lasso, newx = x_test, s = tuning_fit$lambda_1se)),
      RSF = rowSums(forest_pred$chf))
    # Fixed training-derived risk boundaries; evaluated only in the held-out fold.
    train_risk <- as.numeric(predict(lasso, newx = x_train, s = tuning_fit$lambda_1se))
    cuts <- as.numeric(stats::quantile(train_risk, c(1/3, 2/3)))
    risk_group <- if (cuts[1] < cuts[2]) as.character(cut(risks$LASSO_1SE,
      breaks = c(-Inf, cuts, Inf), labels = c("Low", "Medium", "High"))) else rep("Unseparated", nrow(test))
    for (model in names(risks)) {
      i <- length(predictions) + 1L
      predictions[[i]] <- data.frame(repeat_id = r, fold = f, id = test$id,
        time = test$time, status = test$status, death = test$death,
        model, risk = risks[[model]], risk_group)
      fold_metrics[[i]] <- data.frame(repeat_id = r, fold = f, model,
        n_test = nrow(test), deaths = sum(test$death),
        cindex = risk_cindex(test$time, test$death, risks[[model]]))
    }
    co <- as.matrix(coef(lasso, s = tuning_fit$lambda_1se))[, 1]
    selections[[length(selections) + 1L]] <- data.frame(repeat_id = r, fold = f,
      term = names(co), selected = as.integer(co != 0))
    tuning[[length(tuning) + 1L]] <- data.frame(repeat_id = r, fold = f,
      lambda_min = tuning_fit$lambda_min, lambda_1se = tuning_fit$lambda_1se)
  }
}
pred <- do.call(rbind, predictions); metrics <- do.call(rbind, fold_metrics)
save_table(pred, file.path(out, "out_of_fold_predictions.csv"))
save_table(metrics, file.path(out, "fold_metrics.csv"))
save_table(do.call(rbind, tuning), file.path(out, "selected_penalties.csv"))
sel <- stats::aggregate(selected ~ term, do.call(rbind, selections), mean)
save_table(sel, file.path(out, "selection_frequency.csv"))
summary <- do.call(rbind, lapply(split(metrics, metrics$model), function(d)
  data.frame(model = d$model[1], mean_cindex = mean(d$cindex), sd_across_folds = sd(d$cindex),
             min_cindex = min(d$cindex), max_cindex = max(d$cindex), folds = nrow(d))))
save_table(summary, file.path(out, "model_summary.csv"))

# Paired uncertainty across outer test folds. Folds/repeats overlap in training;
# quantiles below describe split sensitivity and are NOT confidence intervals.
wide <- reshape(metrics[c("repeat_id", "fold", "model", "cindex")],
                idvar = c("repeat_id", "fold"), timevar = "model", direction = "wide")
deltas <- do.call(rbind, lapply(c("LASSO", "LASSO_1SE", "RSF"), function(m) {
  delta <- wide[[paste0("cindex.", m)]] - wide$cindex.Cox
  data.frame(model = m, mean_difference_vs_cox = mean(delta),
             split_q025 = unname(quantile(delta, .025)),
             split_q975 = unname(quantile(delta, .975)))
}))
save_table(deltas, file.path(out, "paired_split_differences.csv"))

# Multi-state descriptive estimand: Aalen-Johansen death cumulative incidence.
ci <- cmprsk::cuminc(raw$time, raw$status, cencode = 0)
aj_death <- ci[["1 2"]] # no-group label 1, death status 2
aj_transplant <- ci[["1 1"]]
naive <- survival::survfit(survival::Surv(time, status == 2L) ~ 1, data = raw)
curves <- rbind(
 data.frame(day = aj_death$time, probability = aj_death$est, estimand = "Death cumulative incidence (AJ)"),
 data.frame(day = aj_transplant$time, probability = aj_transplant$est, estimand = "Transplant cumulative incidence (AJ)"),
 data.frame(day = naive$time, probability = 1 - naive$surv, estimand = "1 - KM (transplants censored)"))
save_table(curves, file.path(out, "competing_event_curves.csv"))
five_year <- 5 * 365.25
step_value <- function(t, y, at) if (any(t <= at)) tail(y[t <= at], 1) else 0
save_table(data.frame(day = five_year,
  death_cif = step_value(aj_death$time, aj_death$est, five_year),
  transplant_cif = step_value(aj_transplant$time, aj_transplant$est, five_year),
  naive_one_minus_km = step_value(naive$time, 1-naive$surv, five_year)),
  file.path(out, "five_year_estimands.csv"))

# PH diagnostic is descriptive on the full cohort; it is not an evaluation score.
prep_all <- fit_pbc_preprocessor(raw); all_data <- apply_pbc_preprocessor(raw, prep_all)
full_cox <- survival::coxph(cox_formula, data = all_data, x = TRUE)
ph <- survival::cox.zph(full_cox)
save_table(data.frame(term = rownames(ph$table), ph$table, row.names = NULL), file.path(out, "ph_diagnostic.csv"))

library(ggplot2)
p <- ggplot(metrics, aes(x = model, y = cindex, fill = model)) +
 geom_boxplot(width = .5, alpha = .55, outlier.shape = NA) +
 geom_point(position = position_jitter(width = .07, seed = 4), size = 1.6, alpha = .65) +
 scale_fill_manual(values = c("#17496A", "#228B8D", "#D5913C", "#7C69A4")) +
 coord_cartesian(ylim = c(.5, 1)) + labs(title = "Does model complexity improve discrimination?",
 subtitle = paste(repeats, "repeats x 5 outer folds; LASSO tuned in 4 inner folds; each point is one test fold"),
 x = NULL, y = "Held-out Harrell C (death; transplant censored)") + plot_theme() + theme(legend.position = "none")
save_plot(p, file.path(out, "figures/model_comparison.png"), 9, 5.5)
p <- ggplot(curves, aes(day / 365.25, probability, color = estimand)) + geom_step(linewidth = .9) +
 scale_color_manual(values = c("#D5913C", "#17496A", "#228B8D")) +
 labs(title = "The endpoint changes the probability being estimated",
 subtitle = "Full cohort, descriptive curves; death and transplantation are mutually exclusive first events",
 x = "Years from registration", y = "Probability", color = NULL) + plot_theme()
save_plot(p, file.path(out, "figures/competing_events.png"), 9, 5.5)
p <- ggplot(sel, aes(reorder(term, selected), selected)) + geom_col(fill = "#228B8D", width = .7) +
 coord_flip() + labs(title = "How stable is LASSO feature selection?",
 subtitle = "Fraction of outer fits selecting each term at the training-only one-SE penalty",
 x = NULL, y = "Selection frequency") + plot_theme()
save_plot(p, file.path(out, "figures/selection_stability.png"), 8, 7)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
message("PBC analysis complete. Results: ", out)
