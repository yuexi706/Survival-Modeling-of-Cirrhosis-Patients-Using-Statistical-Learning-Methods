# Behavioral validation. Only checks for projects present in this package are run.
source("R/common.R")
if (file.exists("projects/pbc-survival/R/methods.R")) source("tests/pbc_checks.R")
if (file.exists("projects/bayesian-wine/R/methods.R")) source("tests/wine_checks.R")
message("All requested validation checks passed.")
