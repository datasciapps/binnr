# Input, covariate and outcome preparation ------------------------------------

as_omics <- function(x, call = rlang::caller_env()) {
  if (is.matrix(x) || is.data.frame(x)) {
    x <- list(x = x)
  }
  if (!is.list(x) || length(x) == 0) {
    cli::cli_abort(
      "Omics data must be a matrix, data frame or named list of matrices, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  if (is.null(names(x)) || any(names(x) == "") || anyDuplicated(names(x))) {
    cli::cli_abort("A list of omics matrices must have unique names (e.g. {.val rna}, {.val cnv}).",
                   call = call)
  }
  x <- lapply(x, function(m) {
    if (is.data.frame(m)) m <- as.matrix(m)
    if (!is.numeric(m)) cli::cli_abort("Omics matrices must be numeric.", call = call)
    if (is.null(colnames(m))) {
      cli::cli_abort("Omics matrices need column names (gene identifiers matching the graph).",
                     call = call)
    }
    m
  })
  n <- vapply(x, nrow, integer(1))
  if (length(unique(n)) != 1L) {
    cli::cli_abort(c("All omics matrices must have the same number of rows (samples).",
                     i = "Row counts: {paste(names(n), n, sep = ' = ', collapse = ', ')}."),
                   call = call)
  }
  rn <- lapply(x, rownames)
  if (!any(vapply(rn, is.null, logical(1)))) {
    if (!all(vapply(rn, identical, logical(1), rn[[1]]))) {
      cli::cli_abort("Row names of omics matrices must match (same samples in the same order).",
                     call = call)
    }
  }
  x
}

check_complete <- function(x, what, call = rlang::caller_env()) {
  if (anyNA(x)) {
    cli::cli_abort(c(
      "{what} contain{?s} missing values.",
      i = "Missing data are not supported yet; remove or impute them first."
    ), call = call)
  }
}

# Coerce new data to the omics layers of a fitted network; a single matrix (or
# one-element list) is accepted for single-modality models whatever its name.
as_model_omics <- function(new_data, net, call = rlang::caller_env()) {
  if ((is.matrix(new_data) || is.data.frame(new_data)) && length(net$modalities) == 1L) {
    new_data <- stats::setNames(list(new_data), net$modalities)
  }
  omics <- as_omics(new_data, call = call)
  if (length(omics) == 1L && length(net$modalities) == 1L) names(omics) <- net$modalities
  omics
}

input_matrix <- function(omics, net, scaler = NULL, scale = TRUE,
                         call = rlang::caller_env()) {
  miss <- setdiff(net$modalities, names(omics))
  if (length(miss) > 0) {
    cli::cli_abort("New data lack omics layer{?s} {.val {miss}}.", call = call)
  }
  inputs <- net$nodes[net$nodes$type == "input", ]
  X <- do.call(cbind, lapply(net$modalities, function(m) {
    genes <- inputs$gene[inputs$modality == m]
    mat <- omics[[m]]
    hit <- match(genes, colnames(mat))
    if (anyNA(hit)) {
      cli::cli_abort("Omics layer {.val {m}} is missing {sum(is.na(hit))} feature{?s} used in the model.",
                     call = call)
    }
    mat[, hit, drop = FALSE]
  }))
  storage.mode(X) <- "double"
  check_complete(X, "Omics data", call = call)
  if (is.null(scaler)) {
    if (scale) {
      center <- colMeans(X)
      sds <- sqrt(colSums((X - rep(center, each = nrow(X)))^2) / max(nrow(X) - 1, 1))
      sds[!is.finite(sds) | sds < 1e-8] <- 1
    } else {
      center <- rep(0, ncol(X))
      sds <- rep(1, ncol(X))
    }
    scaler <- list(center = center, scale = sds)
  }
  X <- (X - rep(scaler$center, each = nrow(X))) / rep(scaler$scale, each = nrow(X))
  list(x = X, scaler = scaler)
}

# Covariates enter the prediction head as a scaled design matrix. hardhat
# does the formula/factor bookkeeping (mold() at fit time, forge() at
# prediction time, which checks columns and factor levels); we only add
# centring and scaling, whose parameters are stored in `spec`.
prepare_covariates <- function(data, spec = NULL, scale = TRUE,
                               call = rlang::caller_env()) {
  if (is.null(data) && is.null(spec)) return(list(x = NULL, spec = NULL))
  if (!is.null(spec) && is.null(spec$blueprint)) return(list(x = NULL, spec = spec))
  if (is.null(data)) {
    cli::cli_abort("This model uses covariates; supply them as {.arg covariates}.", call = call)
  }
  if (is.null(spec)) {
    return(mold_covariates(~ ., as.data.frame(data), scale = scale, call = call))
  }
  # forge() errors on missing columns and warns on novel factor levels; both
  # mean the covariates do not match the fit
  forged <- tryCatch(
    hardhat::forge(as.data.frame(data), spec$blueprint),
    error = function(e) {
      cli::cli_abort(c("Covariates are incompatible with those used for fitting.",
                       x = conditionMessage(e)), call = call)
    },
    warning = function(w) {
      cli::cli_abort(c("Covariates are incompatible with those used for fitting.",
                       x = conditionMessage(w)), call = call)
    }
  )
  mm <- design_matrix(forged$predictors)
  check_complete(mm, "Covariates", call = call)
  list(x = scale_columns(mm, spec$center, spec$scale), spec = spec)
}

# Fit-time counterpart: `formula` is a right-hand-side formula
mold_covariates <- function(formula, data, scale = TRUE, call = rlang::caller_env()) {
  rlang::check_installed("hardhat")
  bp <- hardhat::default_formula_blueprint(intercept = TRUE, indicators = "traditional")
  m <- hardhat::mold(formula, data, blueprint = bp)
  mm <- design_matrix(m$predictors)
  if (ncol(mm) == 0) return(list(x = NULL, spec = NULL))
  check_complete(mm, "Covariates", call = call)
  if (scale) {
    center <- colMeans(mm)
    sds <- apply(mm, 2, stats::sd)
    sds[!is.finite(sds) | sds < 1e-8] <- 1
  } else {
    center <- rep(0, ncol(mm))
    sds <- rep(1, ncol(mm))
  }
  spec <- list(blueprint = m$blueprint, columns = colnames(mm), center = center, scale = sds)
  list(x = scale_columns(mm, center, sds), spec = spec)
}

design_matrix <- function(predictors) {
  mm <- as.matrix(predictors)
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

scale_columns <- function(mm, center, scale) {
  (mm - rep(center, each = nrow(mm))) / rep(scale, each = nrow(mm))
}

covariate_names <- function(spec) names(spec$blueprint$ptypes$predictors)

prepare_outcome <- function(y, call = rlang::caller_env()) {
  if (inherits(y, "Surv")) {
    if (attr(y, "type") != "right") {
      cli::cli_abort("Only right-censored {.cls Surv} outcomes are supported.", call = call)
    }
    y <- as.matrix(y)
    check_complete(y, "The outcome", call = call)
    if (sum(y[, "status"]) == 0) {
      cli::cli_abort("The survival outcome has no events.", call = call)
    }
    if (any(y[, "time"] <= 0)) {
      cli::cli_abort("Survival times must be positive.", call = call)
    }
    return(list(family = "cox", time = unname(y[, "time"]),
                status = unname(y[, "status"]), n_out = 1L))
  }
  if (is.character(y) || is.logical(y)) y <- factor(y)
  check_complete(y, "The outcome", call = call)
  if (is.factor(y)) {
    y <- droplevels(y)
    if (nlevels(y) < 2) cli::cli_abort("A factor outcome needs at least two classes.", call = call)
    return(list(family = "classification", y = as.integer(y), levels = levels(y),
                n_out = nlevels(y)))
  }
  if (is.numeric(y)) {
    mu <- mean(y)
    s <- stats::sd(y)
    if (!is.finite(s) || s == 0) s <- 1
    return(list(family = "regression", y = (y - mu) / s, center = mu, scale = s,
                n_out = 1L))
  }
  cli::cli_abort(
    "The outcome must be a factor (classification), numeric (regression) or {.fn survival::Surv} object, not {.obj_type_friendly {y}}.",
    call = call
  )
}

# MultiAssayExperiment support --------------------------------------------------

#' Extract omics matrices from a MultiAssayExperiment
#'
#' Converts selected experiments of a
#' [MultiAssayExperiment::MultiAssayExperiment()] into the sample-by-feature
#' matrices used by [bn_fit()], keeping only patients observed in every
#' selected experiment (one sample per patient), together with the matching
#' clinical data.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param assays Names (or indices) of experiments to use. Defaults to all
#'   experiments that are matrix-like (e.g. `SummarizedExperiment`).
#' @param names Optional short names for the omics layers.
#' @return A list with `omics` (named list of matrices; rows are patients)
#'   and `clinical` (a tibble of `colData`, in the same order).
#' @examplesIf binnr::tcga_cached("BRCA") && requireNamespace("MultiAssayExperiment", quietly = TRUE)
#' brca <- tcga_cohort("BRCA")
#' library(MultiAssayExperiment)
#' mae <- MultiAssayExperiment(
#'   experiments = list(rna = t(brca$omics$rna[1:20, 1:50]),
#'                      cnv = t(brca$omics$cnv[1:20, 1:50])),
#'   colData = S4Vectors::DataFrame(brca$clinical[1:20, ],
#'                                  row.names = brca$clinical$patient[1:20])
#' )
#' d <- mae_to_omics(mae)
#' str(d, max.level = 1)
#' @export
mae_to_omics <- function(mae, assays = NULL, names = NULL) {
  rlang::check_installed(c("MultiAssayExperiment", "SummarizedExperiment"))
  exps <- MultiAssayExperiment::experiments(mae)
  if (is.null(assays)) {
    ok <- vapply(as.list(exps), function(e) {
      is.matrix(e) || methods::is(e, "SummarizedExperiment")
    }, logical(1))
    assays <- base::names(exps)[ok]
  }
  if (is.numeric(assays)) assays <- base::names(exps)[assays]
  sm <- MultiAssayExperiment::sampleMap(mae)
  mats <- lapply(assays, function(a) {
    e <- exps[[a]]
    m <- if (is.matrix(e)) e else as.matrix(SummarizedExperiment::assay(e))
    if (is.null(rownames(m)) && methods::is(e, "SummarizedExperiment")) {
      rd <- SummarizedExperiment::rowData(e)
      sym <- intersect(c("Gene.Symbol", "gene_symbol", "symbol", "gene"), base::names(rd))
      if (length(sym) > 0) rownames(m) <- rd[[sym[[1]]]]
    }
    m <- m[!duplicated(rownames(m)), , drop = FALSE]
    s <- sm[sm$assay == a, ]
    prim <- s$primary[match(colnames(m), s$colname)]
    keep <- !is.na(prim) & !duplicated(prim)
    m <- t(m[, keep, drop = FALSE])
    rownames(m) <- prim[keep]
    storage.mode(m) <- "double"
    m
  })
  base::names(mats) <- names %||% assays
  pts <- Reduce(intersect, lapply(mats, rownames))
  if (length(pts) == 0) cli::cli_abort("No patients are shared by all selected assays.")
  mats <- lapply(mats, function(m) m[pts, , drop = FALSE])
  cd <- as.data.frame(MultiAssayExperiment::colData(mae))
  cd <- cd[pts, , drop = FALSE]
  clinical <- if ("patient" %in% base::names(cd)) {
    tibble::as_tibble(cd)
  } else {
    tibble::as_tibble(cd, rownames = "patient")
  }
  list(omics = mats, clinical = clinical)
}
