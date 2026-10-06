#' @export
print.bn_model <- function(x, ...) {
  fam <- switch(x$family,
    classification = paste0("classification (", length(x$outcome$levels), " classes)"),
    regression = "regression",
    cox = "Cox proportional hazards"
  )
  nodes <- x$network$nodes
  cat_line("<bn_model> biologically informed network for {fam}")
  cat_bullets(c(
    "Omics: {paste0(names(x$network$keep), ' (', lengths(x$network$keep), ')', collapse = ', ')}",
    "Graph: {length(unique(nodes$gene[nodes$type == 'input']))} gene{?s}, {sum(nodes$type == 'pathway')} pathway{?s}, {length(x$network$src)} edge{?s}, depth {max(nodes$level)}",
    "Architecture: {x$arch$preset} preset, embedding dim {x$arch$hidden_dim}, {x$arch$weights} weights, {x$arch$readout} readout",
    "Covariates: {if (x$n_cov > 0) paste(x$n_cov, 'column(s) (late fusion)') else 'none'}",
    "{format(x$n_params, big.mark = ',')} parameters; trained on {x$n_train} of {x$n} samples",
    "Best epoch {x$best_epoch} of {nrow(x$history)} (loss {signif(x$best_loss, 4)}, {round(x$elapsed, 1)}s)"
  ))
  invisible(x)
}

#' @rdname tidy.bn_model
#' @param object A `bn_model`.
#' @details `summary()` returns the same one-row tibble as `glance()`.
#' @export
summary.bn_model <- function(object, ...) {
  glance.bn_model(object)
}

#' Tidy and glance methods for fitted BINNs
#'
#' `tidy()` returns the learned edge weights of the knowledge graph (for
#' `weights = "edge"` architectures) together with node types and levels.
#' `glance()` returns a one-row summary of the fit.
#'
#' Raw edge weights are **not** reliable measures of biological importance
#' (they are not identifiable and vary between runs; see Esser-Skala &
#' Fortelny, 2023). Prefer [node_importance()] and compare across seeds.
#'
#' @param x A `bn_model`.
#' @param ... Unused.
#' @return A tibble.
#' @name tidy.bn_model
#' @examplesIf binnr::tcga_cached("BRCA") && torch::torch_is_installed()
#' brca <- tcga_cohort("BRCA")
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
#' fit <- bn_fit(brca$omics$rna, brca$omics$rna[, "ESR1"], graph = g,
#'               control = bn_control(epochs = 5, seed = 1))
#' tidy(fit)
#' glance(fit)
NULL

#' @rdname tidy.bn_model
#' @exportS3Method generics::tidy
tidy.bn_model <- function(x, ...) {
  net <- x$network
  nodes <- net$nodes
  w <- x$params$edge_weight %||% rep(NA_real_, length(net$src))
  tibble::tibble(
    from = nodes$name[net$src],
    to = nodes$name[net$dst],
    from_type = nodes$type[net$src],
    to_type = nodes$type[net$dst],
    to_level = nodes$level[net$dst],
    to_label = nodes$label[net$dst],
    weight = as.numeric(w)
  )
}

#' @rdname tidy.bn_model
#' @exportS3Method generics::glance
glance.bn_model <- function(x, ...) {
  nodes <- x$network$nodes
  tibble::tibble(
    family = x$family,
    preset = x$arch$preset,
    hidden_dim = x$arch$hidden_dim,
    n = x$n,
    n_train = x$n_train,
    n_features = x$network$n_input,
    n_nodes = nrow(nodes) - x$network$n_input,
    n_edges = length(x$network$src),
    depth = max(nodes$level),
    n_params = x$n_params,
    epochs = nrow(x$history),
    best_epoch = x$best_epoch,
    best_loss = x$best_loss,
    seconds = x$elapsed
  )
}

#' @importFrom generics tidy glance
#' @export
generics::tidy

#' @export
generics::glance

#' Plot training history
#'
#' @param object A `bn_model`.
#' @param ... Unused.
#' @return A ggplot object.
#' @examplesIf binnr::tcga_cached("BRCA") && torch::torch_is_installed() && requireNamespace("ggplot2", quietly = TRUE)
#' brca <- tcga_cohort("BRCA")
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
#' fit <- bn_fit(brca$omics$rna[1:200, ], brca$omics$rna[1:200, "ESR1"], graph = g,
#'               control = bn_control(epochs = 3, seed = 1))
#' ggplot2::autoplot(fit)
#' @exportS3Method ggplot2::autoplot
autoplot.bn_model <- function(object, ...) {
  rlang::check_installed("ggplot2")
  h <- object$history
  long <- rbind(
    data.frame(epoch = h$epoch, set = "training", loss = h$train),
    data.frame(epoch = h$epoch, set = "validation", loss = h$validation)
  )
  long <- long[!is.na(long$loss), ]
  ggplot2::ggplot(long, ggplot2::aes(.data$epoch, .data$loss, colour = .data$set)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = object$best_epoch, linetype = "dashed", colour = "grey50") +
    ggplot2::labs(x = "Epoch", y = "Loss", colour = NULL,
                  title = "Training history",
                  subtitle = paste("Best epoch:", object$best_epoch)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' @export
plot.bn_model <- function(x, ...) {
  print(autoplot.bn_model(x, ...))
  invisible(x)
}
