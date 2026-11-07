# Shared deterministic utilities; run scripts from the repository root.
require_packages <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install missing packages: ", paste(absent, collapse = ", "))
}

stratified_folds <- function(stratum, k, seed) {
  set.seed(seed)
  fold <- integer(length(stratum))
  for (idx in split(seq_along(stratum), stratum)) {
    fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold
}

save_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15),
                   plot.subtitle = ggplot2::element_text(color = "#536477", size = 10),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")
}

save_plot <- function(plot, path, width = 8, height = 5) {
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 180, bg = "white")
}
