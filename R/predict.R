#' Predictions from a fitted BINN
#'
#' Returns predictions as a tibble with one row per sample, using the column
#' naming conventions of tidymodels (`.pred_class`, `.pred_<level>`, `.pred`,
#' `.pred_linear_pred`, `.pred_time`).
#'
#' @param object A `bn_model` from [bn_fit()].
#' @param new_data Omics data for new samples, in the same form as used for
#'   fitting (matrix or named list of matrices with the same feature names).
#'   Defaults to none, which is an error: training data are not stored in the
#'   model.
#' @param covariates For models with covariates: a data frame containing the
#'   covariates (for formula fits, the right-hand-side variables).
#' @param type Type of prediction:
#'   * classification: `"class"` (default) or `"prob"`;
#'   * regression: `"numeric"`;
#'   * survival: `"linear_pred"` (default; the log relative hazard),
#'     `"risk"` (relative hazard) or `"survival"` (survival probabilities at
#'     `eval_time`, from the Breslow baseline hazard, as a list-column).
#' @param eval_time Times at which to predict survival probabilities.
#' @param ... Unused.
#' @return A tibble.
#' @examplesIf binnr::tcga_cached("BRCA") && torch::torch_is_installed()
#' brca <- tcga_cohort("BRCA")
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
#' ok <- which(!is.na(brca$clinical$stage))
#' train <- ok[1:600]
#' test <- ok[-(1:600)]
#' fit <- bn_fit(brca$omics$rna[train, ], brca$clinical$stage[train], graph = g,
#'               control = bn_control(epochs = 5, seed = 1))
#' predict(fit, brca$omics$rna[test, ])
#' predict(fit, brca$omics$rna[test, ], type = "prob")
#' @export
predict.bn_model <- function(object, new_data, covariates = NULL, type = NULL,
                             eval_time = NULL, ...) {
  rlang::check_dots_empty()
  if (missing(new_data)) {
    cli::cli_abort("{.arg new_data} is required; binnr models do not store training data.")
  }
  type <- type %||% switch(object$family,
    classification = "class", regression = "numeric", cox = "linear_pred")
  allowed <- switch(object$family,
    classification = c("class", "prob"),
    regression = "numeric",
    cox = c("linear_pred", "risk", "survival")
  )
  type <- rlang::arg_match0(type, allowed)
  out <- predict_raw(object, new_data, covariates)
  switch(type,
    class = {
      lv <- object$outcome$levels
      tibble::tibble(.pred_class = factor(lv[max.col(out, ties.method = "first")], levels = lv))
    },
    prob = {
      p <- exp(out - apply(out, 1, max))
      p <- p / rowSums(p)
      colnames(p) <- paste0(".pred_", object$outcome$levels)
      tibble::as_tibble(as.data.frame(p, check.names = FALSE))
    },
    numeric = tibble::tibble(.pred = out[, 1] * object$outcome$scale + object$outcome$center),
    linear_pred = tibble::tibble(.pred_linear_pred = out[, 1]),
    risk = tibble::tibble(.pred_risk = exp(out[, 1])),
    survival = {
      if (is.null(eval_time)) cli::cli_abort("{.arg eval_time} is required for {.code type = \"survival\"}.")
      bl <- object$baseline
      H0 <- stats::stepfun(bl$time, c(0, bl$cumhaz))(eval_time)
      rr <- exp(out[, 1] - bl$shift)
      tibble::tibble(.pred = lapply(rr, function(r) {
        tibble::tibble(.eval_time = eval_time, .pred_survival = exp(-H0 * r))
      }))
    }
  )
}

predict_raw <- function(object, new_data, covariates = NULL, module = NULL) {
  omics <- as_model_omics(new_data, object$network)
  X <- input_matrix(omics, object$network, scaler = object$scaler)$x
  C <- new_covariates(object, covariates, nrow(X))
  module <- module %||% rebuild_module(object)
  torch::with_no_grad({
    Xt <- torch::torch_tensor(X, dtype = torch::torch_float())
    Ct <- if (is.null(C)) NULL else torch::torch_tensor(C, dtype = torch::torch_float())
    out <- torch::as_array(module(Xt, Ct))
  })
  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  out
}

new_covariates <- function(object, covariates, n) {
  spec <- object$covariates
  if (is.null(spec) || is.null(spec$blueprint)) {
    if (!is.null(covariates)) {
      cli::cli_warn("This model has no covariates; {.arg covariates} is ignored.")
    }
    return(NULL)
  }
  if (is.null(covariates)) {
    cli::cli_abort("This model uses covariates ({.val {covariate_names(spec)}}); supply them via {.arg covariates}.",
                   call = rlang::caller_env(2))
  }
  C <- prepare_covariates(covariates, spec = spec)$x
  if (nrow(C) != n) {
    cli::cli_abort("{.arg covariates} has {nrow(C)} row{?s} but {.arg new_data} has {n}.",
                   call = rlang::caller_env(2))
  }
  C
}
