#' Node importance scores
#'
#' Scores every node in the network (input features, genes and pathways) by
#' the mean absolute *gradient x activation* of the model output with respect
#' to the node's embedding, a first-order attribution closely related to
#' DeepLIFT and integrated gradients as used with P-NET and BINN.
#'
#' Explanations from BINNs can be unstable across random initialisations
#' (the "Rashomon effect"). Use `seeds` in [bn_importance_stability()] or
#' refit with different seeds before interpreting individual pathways.
#'
#' @param object A `bn_model`.
#' @param new_data Omics data (as for [predict.bn_model()]), typically the
#'   training or test samples.
#' @param covariates Covariates for models that use them.
#' @param target For classification: the class whose logit is explained.
#'   By default, each sample's predicted class. Ignored otherwise.
#' @param normalize If `TRUE`, scores are divided by their maximum within
#'   each node type so that they are comparable across layers.
#' @return A tibble of class `bn_importance` with columns `node`, `label`,
#'   `type`, `level`, `modality`, `importance`, sorted by decreasing
#'   importance.
#' @examplesIf binnr::tcga_cached("BRCA") && torch::torch_is_installed()
#' brca <- tcga_cohort("BRCA")
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
#'                    labels = reactome_tcga$pathways)
#' i <- 1:300
#' fit <- bn_fit(pam50 ~ 1, data = brca$clinical[i, ], omics = list(rna = brca$omics$rna[i, ]),
#'               graph = g, control = bn_control(epochs = 5, seed = 1))
#' imp <- node_importance(fit, list(rna = brca$omics$rna[i, ]))
#' imp[imp$type == "pathway", ]
#' @export
node_importance <- function(object, new_data, covariates = NULL, target = NULL,
                            normalize = FALSE) {
  if (!inherits(object, "bn_model")) cli::cli_abort("{.arg object} must be a {.cls bn_model}.")
  omics <- as_model_omics(new_data, object$network)
  X <- input_matrix(omics, object$network, scaler = object$scaler)$x
  C <- new_covariates(object, covariates, nrow(X))
  module <- rebuild_module(object)
  Xt <- torch::torch_tensor(X, dtype = torch::torch_float(), requires_grad = TRUE)
  Ct <- if (is.null(C)) NULL else torch::torch_tensor(C, dtype = torch::torch_float())
  levels <- module$propagate(Xt)
  for (l in levels[-1]) l$retain_grad()
  levels[[1]]$retain_grad()
  out <- module$readout_head(levels, Ct)
  if (object$family == "classification") {
    if (is.null(target)) {
      cls <- max.col(torch::as_array(out$detach()), ties.method = "first")
    } else {
      cls <- match(target, object$outcome$levels)
      if (is.na(cls)) cli::cli_abort("{.arg target} must be one of {.val {object$outcome$levels}}.")
      cls <- rep(cls, nrow(X))
    }
    sel <- out$gather(2, torch::torch_tensor(matrix(cls, ncol = 1), dtype = torch::torch_long()))
    sel$sum()$backward()
  } else {
    out$sum()$backward()
  }
  score <- unlist(lapply(levels, function(l) {
    ga <- (l$grad * l)$detach()$abs()$sum(dim = 3)$mean(dim = 1)
    as.numeric(torch::as_array(ga))
  }))
  nodes <- object$network$nodes
  res <- tibble::tibble(
    node = nodes$name,
    label = nodes$label,
    type = nodes$type,
    level = nodes$level,
    modality = nodes$modality,
    importance = score
  )
  if (normalize) {
    mx <- stats::ave(res$importance, res$type, FUN = function(v) max(v, 1e-12))
    res$importance <- res$importance / mx
  }
  res <- res[order(-res$importance), ]
  class(res) <- c("bn_importance", class(res))
  res
}

#' Stability of node importance across random seeds
#'
#' Refits a model with several seeds and computes [node_importance()] for
#' each, to assess how reproducible pathway-level explanations are.
#'
#' @param object A `bn_model`, whose original data are re-supplied.
#' @param x,y,covariates The data used to fit `object` (default-method
#'   arguments of [bn_fit()]).
#' @param graph The knowledge graph.
#' @param seeds Integer seeds.
#' @param type Node type to summarise (default `"pathway"`).
#' @return A tibble with one row per node: mean and standard deviation of
#'   importance and mean rank across seeds, plus the per-seed scores in
#'   attribute `"scores"`.
#' @export
bn_importance_stability <- function(object, x, y, graph, covariates = NULL,
                                    seeds = 1:5, type = "pathway") {
  scores <- lapply(seeds, function(s) {
    ctrl <- object$control
    ctrl$seed <- s
    ctrl$verbose <- FALSE
    fit <- suppressMessages(bn_fit(x, y, graph = graph, covariates = covariates,
                                   arch = object$arch, control = ctrl))
    imp <- node_importance(fit, x, covariates = covariates)
    imp <- imp[imp$type == type, ]
    imp$rank <- rank(-imp$importance)
    imp$seed <- s
    imp
  })
  all <- do.call(rbind, lapply(scores, as.data.frame))
  agg <- stats::aggregate(cbind(importance, rank) ~ node + label, data = all,
                          FUN = function(v) c(mean = mean(v), sd = stats::sd(v)))
  res <- tibble::tibble(
    node = agg$node, label = agg$label,
    mean_importance = agg$importance[, "mean"], sd_importance = agg$importance[, "sd"],
    mean_rank = agg$rank[, "mean"], sd_rank = agg$rank[, "sd"]
  )
  res <- res[order(res$mean_rank), ]
  attr(res, "scores") <- tibble::as_tibble(all)
  res
}

#' @rdname node_importance
#' @param n Number of top nodes to show.
#' @param node_type Which node types to show.
#' @param ... Unused.
#' @exportS3Method ggplot2::autoplot
autoplot.bn_importance <- function(object, n = 15, node_type = "pathway", ...) {
  rlang::check_installed("ggplot2")
  d <- object[object$type %in% node_type, ]
  d <- utils::head(d[order(-d$importance), ], n)
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  ggplot2::ggplot(d, ggplot2::aes(.data$importance, .data$label)) +
    ggplot2::geom_col(fill = "#3B6EA8") +
    ggplot2::labs(x = "Mean |gradient x activation|", y = NULL,
                  title = paste("Top", nrow(d), paste(node_type, collapse = "/"), "nodes")) +
    ggplot2::theme_minimal()
}

#' Node activations (pathway activity scores)
#'
#' Returns the hidden representation of every gene and pathway node for each
#' sample. With the default BINN architecture (`hidden_dim = 1`) this is one
#' interpretable "activity score" per pathway and sample, which can be passed
#' to standard statistical tools: differential activity tests with
#' `stats::kruskal.test()` or limma, survival analysis with
#' `survival::coxph()`, heatmaps, and so on.
#'
#' Activation scores are defined only up to sign and scale (a pathway's
#' downstream weights can flip or rescale it), so compare them within a fitted
#' model and check stability across seeds.
#'
#' @inheritParams node_importance
#' @param type Node types to return (default genes and pathways).
#' @param format `"wide"` (a samples x nodes matrix; for `hidden_dim > 1`
#'   columns are `node[k]`) or `"long"` (a tibble with one row per sample, node
#'   and dimension).
#' @return A matrix or tibble; see `format`.
#' @examplesIf binnr::tcga_cached("BRCA") && torch::torch_is_installed()
#' brca <- tcga_cohort("BRCA")
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
#'                    labels = reactome_tcga$pathways)
#' i <- 1:300
#' fit <- bn_fit(pam50 ~ 1, data = brca$clinical[i, ], omics = list(rna = brca$omics$rna[i, ]),
#'               graph = g, control = bn_control(epochs = 5, seed = 1))
#' act <- node_activations(fit, list(rna = brca$omics$rna[i, ]), type = "pathway")
#' dim(act)
#' @export
node_activations <- function(object, new_data, covariates = NULL,
                             type = c("gene", "pathway"),
                             format = c("wide", "long")) {
  if (!inherits(object, "bn_model")) cli::cli_abort("{.arg object} must be a {.cls bn_model}.")
  format <- rlang::arg_match(format)
  omics <- as_model_omics(new_data, object$network)
  X <- input_matrix(omics, object$network, scaler = object$scaler)$x
  module <- rebuild_module(object)
  H <- torch::with_no_grad({
    levels <- module$propagate(torch::torch_tensor(X, dtype = torch::torch_float()))
    torch::as_array(torch::torch_cat(levels, dim = 2))
  })
  nodes <- object$network$nodes
  keep <- which(nodes$type %in% type)
  d <- dim(H)[3]
  sample_ids <- rownames(omics[[1]]) %||% as.character(seq_len(nrow(X)))
  if (format == "wide") {
    M <- matrix(H[, keep, , drop = FALSE], nrow = nrow(X))
    cn <- if (d == 1) nodes$name[keep] else
      paste0(rep(nodes$name[keep], times = d), "[", rep(seq_len(d), each = length(keep)), "]")
    dimnames(M) <- list(sample_ids, cn)
    attr(M, "labels") <- stats::setNames(nodes$label[keep], nodes$name[keep])
    return(M)
  }
  tibble::tibble(
    sample = rep(sample_ids, times = length(keep) * d),
    node = rep(rep(nodes$name[keep], each = nrow(X)), times = d),
    label = rep(rep(nodes$label[keep], each = nrow(X)), times = d),
    type = rep(rep(nodes$type[keep], each = nrow(X)), times = d),
    dim = rep(seq_len(d), each = nrow(X) * length(keep)),
    activation = as.numeric(H[, keep, , drop = FALSE])
  )
}
