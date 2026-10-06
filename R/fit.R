#' Fit a biologically informed neural network
#'
#' Fits a message-passing neural network whose computational graph is a
#' pathway knowledge graph. Omics features (e.g. RNA expression, copy number,
#' mutations) are mapped to genes, genes to pathways and pathways to parent
#' pathways; predictions are made from pathway embeddings, optionally
#' together with clinical covariates (late fusion).
#'
#' The type of model is determined by the outcome:
#'
#' * a **factor** (or character/logical) gives classification with a softmax
#'   head and cross-entropy loss;
#' * a **numeric** vector gives regression (squared error);
#' * a [survival::Surv()] object gives a Cox proportional hazards model
#'   trained by maximising the partial likelihood.
#'
#' The **formula method** is the most convenient interface: the left-hand side
#' is the outcome, the right-hand side lists clinical covariates (`~ 1` for
#' none), `data` holds these variables and `omics` the omics layers, with rows
#' in the same order. Rows with missing values are dropped with a message. If
#' `data` is a `MultiAssayExperiment`, omics layers and `colData` are
#' extracted automatically (see [mae_to_omics()]).
#'
#' @param x A formula, a numeric matrix / data frame of omics features
#'   (samples in rows, genes in columns), or a named list of such matrices
#'   (one per omics layer).
#' @param y Outcome for the default method: factor, numeric or `Surv`.
#' @param graph A [bn_graph] (e.g. from [pathway_graph()] or
#'   [reactome_graph()]). It is pruned to the observed genes.
#' @param covariates Optional data frame of non-omics covariates (default
#'   method).
#' @param data A data frame (or `MultiAssayExperiment`) containing the
#'   variables in the formula.
#' @param omics Omics matrix or named list of matrices (formula method).
#' @param assays For `MultiAssayExperiment` data: which experiments to use.
#' @param arch Architecture, see [bn_arch()].
#' @param control Training options, see [bn_control()].
#' @param ... Passed between methods.
#'
#' @return An object of class `bn_model`, with methods for [predict()],
#'   [print()], [summary()], [tidy()][tidy.bn_model()], [glance()][glance.bn_model()],
#'   [autoplot()][autoplot.bn_model()] and [node_importance()].
#'   Fitted parameters are stored as plain R arrays, so models can be saved
#'   with [saveRDS()].
#' @seealso [bn_arch()], [bn_control()], [predict.bn_model()],
#'   [node_importance()]
#' @examplesIf binnr::tcga_cached("BRCA") && torch::torch_is_installed()
#' brca <- tcga_cohort("BRCA")
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
#'                    labels = reactome_tcga$pathways)
#'
#' # PAM50 subtype from RNA-seq and copy number, adjusting for age
#' # (a few epochs on a subset, to keep the example fast)
#' i <- 1:300
#' fit <- bn_fit(pam50 ~ age, data = brca$clinical[i, ],
#'               omics = lapply(brca$omics[c("rna", "cnv")], function(m) m[i, ]),
#'               graph = g, control = bn_control(epochs = 5, seed = 1))
#' fit
#' \donttest{
#' # Overall survival, relaxed GNN architecture, on the full cohort
#' library(survival)
#' fit_os <- bn_fit(Surv(os_time, os_event) ~ age + stage, data = brca$clinical,
#'                  omics = brca$omics, graph = g, arch = bn_arch("gnn"),
#'                  control = bn_control(epochs = 20, seed = 1))
#' fit_os
#' }
#' @export
bn_fit <- function(x, ...) {
  UseMethod("bn_fit")
}

#' @rdname bn_fit
#' @export
bn_fit.default <- function(x, y, graph, covariates = NULL,
                           arch = bn_arch(), control = bn_control(), ...) {
  rlang::check_dots_empty()
  if (inherits(x, "MultiAssayExperiment")) {
    cli::cli_abort("For a {.cls MultiAssayExperiment}, use the formula interface: {.code bn_fit(outcome ~ covariates, data = mae, graph = g)}.")
  }
  check_graph(graph)
  if (!inherits(arch, "bn_arch")) cli::cli_abort("{.arg arch} must be created by {.fn bn_arch}.")
  if (!inherits(control, "bn_control")) cli::cli_abort("{.arg control} must be created by {.fn bn_control}.")
  omics <- as_omics(x)
  n <- nrow(omics[[1]])
  yi <- prepare_outcome(y)
  n_y <- if (yi$family == "cox") length(yi$time) else length(yi$y)
  if (n_y != n) cli::cli_abort("{.arg y} has {n_y} value{?s} but the omics data have {n} row{?s}.")
  cov <- prepare_covariates(covariates, scale = control$scale)
  if (!is.null(cov$x) && nrow(cov$x) != n) {
    cli::cli_abort("{.arg covariates} has {nrow(cov$x)} row{?s} but the omics data have {n}.")
  }
  fit_bn_impl(omics, yi, cov, graph, arch, control, call = generic_call(match.call()))
}

#' @rdname bn_fit
#' @export
bn_fit.formula <- function(x, data, omics = NULL, graph, assays = NULL,
                           arch = bn_arch(), control = bn_control(), ...) {
  rlang::check_dots_empty()
  formula <- x
  if (inherits(data, "MultiAssayExperiment")) {
    d <- mae_to_omics(data, assays = assays)
    if (is.null(omics)) omics <- d$omics
    data <- d$clinical
  }
  if (is.null(omics)) cli::cli_abort("{.arg omics} must be supplied.")
  omics <- as_omics(omics)
  data <- as.data.frame(data)
  if (nrow(data) != nrow(omics[[1]])) {
    cli::cli_abort("{.arg data} has {nrow(data)} row{?s} but {.arg omics} has {nrow(omics[[1]])}; rows must correspond to the same samples.")
  }
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  tt <- stats::delete.response(stats::terms(mf))
  rhs <- mf[, -1, drop = FALSE]
  check_graph(graph)
  graph_genes <- bn_nodes(graph)$name[bn_nodes(graph)$type == "gene"]
  ok <- stats::complete.cases(rhs) & !is_na_outcome(y) &
    Reduce(`&`, lapply(omics, function(m) {
      stats::complete.cases(m[, colnames(m) %in% graph_genes, drop = FALSE])
    }))
  if (sum(ok) < 10L) {
    cli::cli_abort("Only {sum(ok)} sample{?s} ha{?s/ve} complete data; at least 10 are needed.")
  }
  if (!all(ok)) {
    cli::cli_inform("Dropping {sum(!ok)} sample{?s} with missing values, leaving {sum(ok)}.")
    omics <- lapply(omics, function(m) m[ok, , drop = FALSE])
    y <- subset_outcome(y, ok)
    data <- data[ok, , drop = FALSE]
  }
  yi <- prepare_outcome(y)
  cov <- if (length(attr(tt, "term.labels")) > 0) {
    mold_covariates(stats::formula(tt), data, scale = control$scale)
  } else {
    list(x = NULL, spec = NULL)
  }
  fit_bn_impl(omics, yi, cov, graph, arch, control, call = generic_call(match.call()),
              formula = formula)
}

# Record calls as `bn_fit(...)` so that update() works on exported names
generic_call <- function(call) {
  call[[1]] <- quote(bn_fit)
  call
}

is_na_outcome <- function(y) {
  if (inherits(y, "Surv")) return(rowSums(is.na(as.matrix(y))) > 0)
  is.na(y)
}

subset_outcome <- function(y, i) {
  if (inherits(y, "Surv")) return(y[i])
  y[i]
}

# Core implementation -----------------------------------------------------------

fit_bn_impl <- function(omics, yi, cov, graph, arch, control, call, formula = NULL) {
  rlang::check_installed("torch")
  if (!torch::torch_is_installed()) {
    cli::cli_abort(c("The torch backend (libtorch) is not installed.",
                     i = "Run {.run torch::install_torch()} once."))
  }
  gene_layer <- arch$gene_layer %||% (length(omics) > 1)
  net <- compile_network(graph, lapply(omics, colnames), gene_layer)
  n_used <- lengths(net$keep)
  n_all <- vapply(omics, ncol, integer(1))
  if (sum(n_used) == 0) {
    cli::cli_abort("None of the omics features match genes in {.arg graph}.")
  }
  if (control$verbose || any(n_used < n_all)) {
    msg <- c("i" = "Using {sum(n_used)} of {sum(n_all)} omics feature{?s} that map to the knowledge graph.")
    if (length(omics) > 1) {
      msg <- c(msg, " " = paste(names(n_used), n_used, sep = ": ", collapse = ", "))
    }
    cli::cli_inform(msg)
  }
  inp <- input_matrix(omics, net, scale = control$scale)
  X <- inp$x
  C <- cov$x
  n <- nrow(X)

  seed <- control$seed
  if (!is.null(seed)) {
    withr::local_seed(seed)
    local_torch_seed(seed)
  }

  # Train / validation split
  val_idx <- integer()
  if (control$validation > 0) {
    strata <- if (yi$family == "classification") yi$y else if (yi$family == "cox") yi$status else rep(1L, n)
    val_idx <- unlist(lapply(split(seq_len(n), strata), function(i) {
      k <- round(length(i) * control$validation)
      i[sample.int(length(i), k)]
    }), use.names = FALSE)
  }
  train_idx <- setdiff(seq_len(n), val_idx)

  device <- torch::torch_device(control$device)
  flt <- function(m) torch::torch_tensor(m, dtype = torch::torch_float(), device = device)
  Xt <- flt(X)
  Ct <- if (is.null(C)) NULL else flt(C)
  target <- make_target(yi, device)

  module <- bn_module(net, arch, n_out = yi$n_out, n_cov = ncol(C) %||% 0L,
                      out_bias = yi$family != "cox")
  module$to(device = device)
  opt <- torch::optim_adamw(module$parameters, lr = control$lr,
                            weight_decay = control$weight_decay)

  batch_size <- control$batch_size %||% if (yi$family == "cox") length(train_idx) else 64L
  batch_size <- min(batch_size, length(train_idx))
  loss_fn <- loss_function(yi$family)

  eval_loss <- function(idx) {
    module$eval()
    torch::with_no_grad({
      out <- module(Xt[idx, , drop = FALSE], if (is.null(Ct)) NULL else Ct[idx, , drop = FALSE])
      loss_fn(out, target, idx)$item()
    })
  }

  history <- vector("list", control$epochs)
  best <- list(loss = Inf, epoch = 0L, state = NULL)
  wait <- 0L
  t0 <- Sys.time()
  for (epoch in seq_len(control$epochs)) {
    module$train()
    perm <- train_idx[sample.int(length(train_idx))]
    batches <- split(perm, ceiling(seq_along(perm) / batch_size))
    nb <- length(batches)
    if (nb > 1L && length(batches[[nb]]) < 2L) {
      batches[[nb - 1L]] <- c(batches[[nb - 1L]], batches[[nb]])
      batches[[nb]] <- NULL
    }
    tr_loss <- 0
    for (b in batches) {
      if (length(b) < 2L) next
      opt$zero_grad()
      out <- module(Xt[b, , drop = FALSE], if (is.null(Ct)) NULL else Ct[b, , drop = FALSE])
      loss <- loss_fn(out, target, b)
      total <- loss
      if (control$l1 > 0 && arch$weights == "edge") {
        total <- total + control$l1 * module$edge_weight$abs()$sum()
      }
      total$backward()
      opt$step()
      tr_loss <- tr_loss + loss$item() * length(b)
    }
    tr_loss <- tr_loss / length(perm)
    # R's garbage collector does not see torch's memory; free tensors promptly
    gc(verbose = FALSE, full = FALSE)
    va_loss <- if (length(val_idx) > 1) eval_loss(val_idx) else NA_real_
    history[[epoch]] <- c(epoch = epoch, train = tr_loss, validation = va_loss)
    monitor <- if (is.na(va_loss)) tr_loss else va_loss
    if (is.finite(monitor) && monitor < best$loss - 1e-6) {
      best <- list(loss = monitor, epoch = epoch, state = get_params(module))
      wait <- 0L
    } else {
      wait <- wait + 1L
    }
    if (control$verbose && (epoch %% 10 == 0 || epoch == 1)) {
      cli::cli_inform("Epoch {epoch}: train {signif(tr_loss, 4)}, validation {signif(va_loss, 4)}")
    }
    if (!is.finite(tr_loss)) {
      cli::cli_warn("Training diverged at epoch {epoch}; try a smaller learning rate.")
      break
    }
    if (wait >= control$patience) break
  }
  history <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  history <- tibble::as_tibble(as.data.frame(history))
  if (is.null(best$state)) best$state <- get_params(module)
  set_params(module, best$state)

  object <- structure(list(
    call = call,
    formula = formula,
    family = yi$family,
    outcome = yi[setdiff(names(yi), c("y", "time", "status"))],
    arch = arch,
    control = control,
    network = net,
    params = best$state,
    scaler = inp$scaler,
    covariates = cov$spec,
    n_cov = ncol(C) %||% 0L,
    history = history,
    best_epoch = best$epoch,
    best_loss = best$loss,
    n = n,
    n_train = length(train_idx),
    elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    n_params = sum(vapply(best$state, length, integer(1)))
  ), class = "bn_model")

  if (yi$family == "cox") {
    module$eval()
    lp <- torch::with_no_grad(as.numeric(torch::as_array(module(Xt, Ct)$cpu())))
    object$baseline <- breslow(lp, yi$time, yi$status)
  }
  object
}

make_target <- function(yi, device) {
  switch(yi$family,
    classification = list(y = torch::torch_tensor(yi$y, dtype = torch::torch_long(), device = device)),
    regression = list(y = torch::torch_tensor(yi$y, dtype = torch::torch_float(), device = device)),
    cox = list(time = yi$time,
               status = torch::torch_tensor(yi$status, dtype = torch::torch_float(), device = device))
  )
}

loss_function <- function(family) {
  switch(family,
    classification = function(out, target, idx) {
      torch::nnf_cross_entropy(out, target$y[idx])
    },
    regression = function(out, target, idx) {
      torch::nnf_mse_loss(out$squeeze(2), target$y[idx])
    },
    cox = function(out, target, idx) {
      cox_partial_loss(out$squeeze(2), target$time[idx], target$status[idx])
    }
  )
}

# Negative Cox partial log-likelihood (Breslow approximation for ties)
cox_partial_loss <- function(eta, time, status) {
  ord <- order(time, decreasing = TRUE)
  ts <- time[ord]
  # With times sorted in decreasing order, the risk set of a subject is
  # everyone up to the *last* subject tied with it.
  last <- stats::ave(seq_along(ts), ts, FUN = max)
  dev <- eta$device
  ord_t <- torch::torch_tensor(ord, dtype = torch::torch_long(), device = dev)
  eta <- eta$index_select(1, ord_t)
  status <- status$index_select(1, ord_t)
  log_risk <- torch::torch_logcumsumexp(eta, dim = 1)$index_select(
    1, torch::torch_tensor(last, dtype = torch::torch_long(), device = dev)
  )
  n_events <- torch::torch_clamp(status$sum(), min = 1)
  -((eta - log_risk) * status)$sum() / n_events
}

# Breslow baseline cumulative hazard, via survival::basehaz() on a Cox model
# with the fitted linear predictor as an offset (no free coefficients).
breslow <- function(lp, time, status) {
  rlang::check_installed("survival")
  fit <- survival::coxph(survival::Surv(time, status) ~ offset(lp), ties = "breslow")
  bh <- survival::basehaz(fit, centered = FALSE)
  list(time = bh$time, cumhaz = bh$hazard, shift = 0)
}

get_params <- function(module) {
  p <- module$parameters
  lapply(p, function(t) {
    a <- torch::as_array(t$detach()$cpu())
    if (is.null(dim(a))) a <- array(a, dim = length(a))
    a
  })
}

set_params <- function(module, params) {
  p <- module$parameters
  torch::with_no_grad({
    for (nm in names(params)) {
      p[[nm]]$copy_(torch::torch_tensor(params[[nm]], dtype = p[[nm]]$dtype))
    }
  })
  invisible(module)
}

# Set torch's RNG seed for the calling frame, restoring the previous state
local_torch_seed <- function(seed = NULL, .local_envir = parent.frame()) {
  state <- torch::torch_get_rng_state()
  withr::defer(torch::torch_set_rng_state(state), envir = .local_envir)
  if (!is.null(seed)) torch::torch_manual_seed(seed)
  invisible(state)
}

# Rebuild a torch module from a fitted model (without disturbing either RNG)
rebuild_module <- function(object, device = "cpu") {
  local_torch_seed()
  module <- withr::with_preserve_seed(
    bn_module(object$network, object$arch, n_out = object$outcome$n_out,
              n_cov = object$n_cov, out_bias = object$family != "cox")
  )
  set_params(module, object$params)
  module$to(device = device)
  module$eval()
  module
}
