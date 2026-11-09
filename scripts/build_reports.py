"""Generate readable result summaries from committed analysis CSVs; standard library only."""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def table(records, fields):
    lines = ["| " + " | ".join(label for _, label in fields) + " |",
             "| " + " | ".join("---" for _ in fields) + " |"]
    for record in records:
        cells = []
        for key, _ in fields:
            value = record[key]
            try:
                value = f"{float(value):.4f}"
            except ValueError:
                pass
            cells.append(value)
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


pbc = ROOT / "projects/pbc-survival"
if pbc.exists():
    result = pbc / "results"
    summary = rows(result / "model_summary.csv")
    cif = rows(result / "five_year_estimands.csv")[0]
    ph = {d["term"]: d for d in rows(result / "ph_diagnostic.csv")}
    content = """# PBC extension results

These are executed extension results from `survival::pbc`, not a reproduction of the missing course CSV. The cohort has 418 records: 232 administratively censored, 25 transplants and 161 deaths.

## Repeated outer-fold discrimination

""" + table(summary, [("model", "Model"), ("mean_cindex", "Mean C"),
                           ("sd_across_folds", "Fold SD"), ("min_cindex", "Minimum"),
                           ("max_cindex", "Maximum")]) + """

Each value summarizes 15 test folds from three repetitions of five-fold evaluation. The SD describes split variation, not an inferential standard error. LASSO uses four inner folds with imputation refitted in each training fold. The Cox model has six prespecified clinical predictors; LASSO and the fixed-parameter forest use the broader 16-feature set.

LASSO has the largest mean C, with an average paired advantage over Cox of about 0.0082. Its paired split differences vary in sign. The one-SE LASSO is close to Cox, and the forest does not improve the mean in this experiment. This is evidence of modest, split-sensitive differences, not established algorithmic superiority or clinical benefit.

## Competing-event probability

""" + f"""At five years, full-cohort Aalen-Johansen death cumulative incidence is **{float(cif['death_cif']):.2%}**, transplant cumulative incidence is **{float(cif['transplant_cif']):.2%}**, and `1 - KM` with transplants censored is **{float(cif['naive_one_minus_km']):.2%}**. The estimands differ even when the numerical gap is small. These curves are descriptive; they do not validate calibrated individual death predictions.

## Model checking

The full-cohort descriptive Schoenfeld-residual diagnostic has global p = **{float(ph['GLOBAL']['p']):.4f}** and prothrombin-time p = **{float(ph['protime']['p']):.4f}**. This suggests investigating time-varying effects rather than assuming the Cox structure is adequate. These are exploratory, unadjusted diagnostic tests after single imputation, not a final inferential analysis.

## Defensible conclusion

The most useful extension is a validation design that separates model development from evaluation and makes the prediction target explicit. Small apparent gains need to be assessed against split sensitivity, missing-data assumptions and model diagnostics. No external validation, multiple imputation or time-varying model was completed.

![Held-out evaluation](results/figures/model_comparison.png)

![Competing events](results/figures/competing_events.png)
"""
    (pbc / "RESULTS.md").write_text(content)

wine = ROOT / "projects/bayesian-wine"
if wine.exists():
    result = wine / "results"
    summary = rows(result / "model_summary.csv")
    paired = rows(result / "paired_group_bootstrap.csv")
    diagnostics = rows(result / "mcmc_diagnostics.csv")
    sampler = rows(result / "sampler_summary.csv")
    max_rhat = max(float(d["rhat"]) for d in diagnostics if d["rhat"])
    bulk = min(float(d["ess_bulk"]) for d in diagnostics if d["ess_bulk"])
    tail = min(float(d["ess_tail"]) for d in diagnostics if d["ess_tail"])
    ebfmi = min(float(d["ebfmi"]) for d in sampler)
    content = """# Bayesian wine extension results

The UCI red-wine extension uses 1,599 records: 1,285 training rows in 1,089 feature groups and 314 held-out rows in 270 groups. Identical predictor vectors do not cross partitions. These are new extension results, not the original 1,143-row course experiment.

## Held-out predictive scores

Lower values are better. The Gaussian distribution is discretized onto the same six rating categories before scoring. Absolute error refers to the expected score under these category probabilities.

""" + table(summary, [("model", "Model"), ("negative_log_score", "Negative log score"),
                           ("rps", "Normalized RPS"), ("absolute_error", "Absolute error")]) + """

The primary ordinal model improves mean log score and RPS modestly, while its expected-rating absolute error is slightly worse. The sensitivity prior produces similar scores; it is not selected using the held-out results. A more suitable response likelihood need not win on every predictive criterion.

## Paired uncertainty

""" + table(paired, [("metric", "Metric"),
                        ("difference_ordinal_minus_gaussian", "Ordinal minus Gaussian"),
                        ("conditional_ci_lower", "Conditional 95% lower"),
                        ("conditional_ci_upper", "Conditional 95% upper")]) + """

Intervals use 1,000 bootstrap draws of held-out feature groups, with model fits and partition fixed. They do not include refitting uncertainty. Small absolute differences on one split are not proof of broad superiority.

## Computation and predictive checking

""" + f"""Across the three posterior fits, maximum rank-normalized split R-hat is **{max_rhat:.4f}**, minimum bulk ESS is **{bulk:.0f}**, and minimum tail ESS is **{tail:.0f}**. There are **{sum(int(d['divergences']) for d in sampler)}** post-warmup divergences, maximum tree depth is **{max(int(d['maximum_treedepth']) for d in sampler)}**, and minimum chain E-BFMI is **{ebfmi:.3f}**. The separate prior-predictive sampler also passes the recorded checks.

Ordinal posterior predictive frequencies better reproduce the concentrated observed rating distribution than the discretized Gaussian baseline in the training check. This does not by itself validate the proportional-odds assumption. Prior sensitivity intervals and a prior predictive frequency check are saved alongside the held-out scores.

## Defensible conclusion

The project illustrates how likelihood choice changes predictive distributions and the conclusions supported by different loss functions. It supports careful checking and reporting, not a new algorithm, a causal claim about chemistry, or a universally superior predictor.

![Held-out metric comparison](results/figures/heldout_comparison.png)

![Prior sensitivity](results/figures/prior_sensitivity.png)
"""
    (wine / "RESULTS.md").write_text(content)

print("Result summaries generated from analysis CSVs.")
