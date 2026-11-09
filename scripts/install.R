options(repos=c(CRAN="https://cloud.r-project.org"))
packages <- c("survival", "glmnet", "ranger", "cmprsk", "ggplot2", "renv")
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) install.packages(missing)
message("Dependencies available; use renv::restore() for locked versions.")
