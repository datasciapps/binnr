skip_if_no_torch <- function() {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch not installed")
}

# A small layered hierarchy: 4 genes -> 2 pathways -> 1 root
toy_layered <- function() {
  bn_graph(data.frame(
    from = c("G1", "G2", "G2", "G3", "G4", "P1", "P2"),
    to   = c("P1", "P1", "P2", "P2", "P2", "P3", "P3")
  ))
}

# A ragged DAG with a skip edge (G3 -> P3) and two roots
toy_ragged <- function() {
  bn_graph(data.frame(
    from = c("G1", "G2", "G2", "G3", "P1", "P2", "G4", "P2"),
    to   = c("P1", "P1", "P2", "P3", "P3", "P3", "P4", "P4")
  ))
}

toy_data <- function(n = 40, genes = paste0("G", 1:4), seed = 1) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * length(genes)), n, dimnames = list(NULL, genes))
    y <- x[, 1] - x[, 3] + stats::rnorm(n, sd = 0.1)
  })
  list(x = x, y = y)
}

quiet_fit <- function(...) suppressMessages(bn_fit(...))

# TCGA example cohort: only available once cached (see ?tcga_cohort)
brca_or_skip <- function() {
  testthat::skip_if_not(tcga_cached("BRCA"), "TCGA BRCA cohort not cached; run tcga_cohort(\"BRCA\")")
  tcga_cohort("BRCA")
}
