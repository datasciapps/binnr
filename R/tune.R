#' Choose an architecture and graph depth by validation loss
#'
#' Fits [bn_fit()] once per row of a small grid and keeps the fit with the
#' lowest validation loss (the loss used for early stopping, computed on the
#' `validation` fraction held out from the training data; see
#' [bn_control()]). Nothing outside the training data is used, so the result
#' can be evaluated on a separate test set or inside cross-validation without
#' leakage.
#'
#' The grid can vary anything [bn_arch()] accepts (`hidden_dim`, `weights`,
#' `readout`, `dropout`, ...), the training options `lr`, `weight_decay` and
#' `l1`, and `max_level`, the number of pathway levels kept from the knowledge
#' graph (see [bn_prune()]). Feedforward BINNs such as P-NET fix this depth
#' (five levels); binnr uses the whole hierarchy by default (`Inf`).
#'
#' @inheritParams bn_fit
#' @param grid A data frame with one row per configuration, e.g. from
#'   [expand.grid()]. Columns are matched by name to arguments of [bn_arch()],
#'   to `lr`, `weight_decay`, `l1` of [bn_control()], and to `max_level`.
#'   Missing columns take their values from `arch` and `control`.
#' @param x,y,covariates As for the default method of [bn_fit()].
#' @param verbose Print the validation loss of each configuration.
#' @return The best `bn_model`, with an extra element `tuning`: a tibble of the
#'   grid with the validation loss, training epochs and seconds of each fit.
#' @examplesIf binnr::tcga_cached("BRCA") && torch::torch_is_installed()
#' brca <- tcga_cohort("BRCA")
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
#' i <- which(!is.na(brca$clinical$pam50))[1:200]
#' fit <- bn_tune(brca$omics$rna[i, ], brca$clinical$pam50[i], graph = g,
#'                grid = expand.grid(hidden_dim = c(1, 4), max_level = c(5, Inf)),
#'                control = bn_control(epochs = 3, seed = 1))
#' fit$tuning
#' @export
bn_tune <- function(x, y, graph, grid, covariates = NULL,
                    arch = bn_arch(), control = bn_control(), verbose = FALSE) {
  check_graph(graph)
  if (!is.data.frame(grid) || nrow(grid) == 0) {
    cli::cli_abort("{.arg grid} must be a data frame with at least one row.")
  }
  arch_args <- c("hidden_dim", "weights", "readout", "activation", "dropout", "gene_layer")
  ctrl_args <- c("lr", "weight_decay", "l1")
  unknown <- setdiff(names(grid), c(arch_args, ctrl_args, "max_level", "preset"))
  if (length(unknown)) {
    cli::cli_abort("Unknown grid column{?s} {.val {unknown}}.")
  }
  fits <- vector("list", nrow(grid))
  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    row <- as.list(grid[i, , drop = FALSE])
    a <- arch
    for (nm in intersect(names(row), arch_args)) a[[nm]] <- row[[nm]]
    if (!is.null(row$preset)) {
      a <- bn_arch(as.character(row$preset))
      for (nm in intersect(names(row), arch_args)) a[[nm]] <- row[[nm]]
    }
    a$hidden_dim <- as.integer(a$hidden_dim)
    ctrl <- control
    for (nm in intersect(names(row), ctrl_args)) ctrl[[nm]] <- row[[nm]]
    ctrl$verbose <- FALSE
    g <- if (!is.null(row$max_level) && is.finite(row$max_level)) {
      bn_prune(graph, max_level = row$max_level)
    } else {
      graph
    }
    t0 <- Sys.time()
    fit <- suppressMessages(bn_fit(x, y, graph = g, covariates = covariates,
                                   arch = a, control = ctrl))
    fits[[i]] <- fit
    results[[i]] <- tibble::tibble(
      config = i, validation_loss = fit$best_loss, epochs = nrow(fit$history),
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    if (verbose) {
      cli::cli_inform("Config {i}/{nrow(grid)}: validation loss {signif(fit$best_loss, 4)}")
    }
  }
  tuning <- dplyr_free_bind(grid, do.call(rbind, results))
  best <- which.min(tuning$validation_loss)
  out <- fits[[best]]
  out$tuning <- tuning
  out$tuning$selected <- seq_len(nrow(tuning)) == best
  out
}

dplyr_free_bind <- function(grid, res) {
  g <- tibble::as_tibble(grid)
  g$config <- seq_len(nrow(g))
  tibble::as_tibble(merge(g, res, by = "config"))
}
