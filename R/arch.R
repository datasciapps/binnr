#' Network architecture specification
#'
#' Describes how information flows through the knowledge graph. The defaults
#' reproduce a classic feedforward BINN (one scalar node per pathway, a
#' distinct weight on every edge, predictions from the top-level pathways
#' only); each argument relaxes one of the constraints discussed in
#' "Are BINNs just GNNs with extra steps?".
#'
#' Message passing follows a *topological* schedule: each pathway node is
#' updated once, after all of its children, as
#' \deqn{h_v = \sigma\left(\sum_{u \to v} w_{uv}\, h_u W_{\ell(v)} + b_v\right),}
#' where `hidden_dim` sets the dimension of \eqn{h_v}, `weights = "edge"`
#' learns a separate scalar \eqn{w_{uv}} for each edge while
#' `weights = "shared"` fixes \eqn{w_{uv} = 1/\mathrm{indeg}(v)} (mean
#' aggregation, as in a GCN), and \eqn{W_\ell} is a `hidden_dim` x
#' `hidden_dim` matrix shared by all nodes at level \eqn{\ell} (omitted when
#' `hidden_dim = 1` and `weights = "edge"`, which is exactly a masked sparse
#' linear layer). Because updates follow the DAG, edges that skip levels need
#' no padding nodes. To reproduce the *stratified* schedule of feedforward
#' BINNs (P-NET, BINN), which updates one padded stratum at a time, fit a
#' graph prepared with [bn_pad()] instead.
#'
#' @param preset Optional starting point: `"binn"` (the defaults) or `"gnn"`
#'   (4-dimensional node embeddings with GCN-style shared weights, i.e. mean
#'   aggregation followed by a learnable transform per level). Other
#'   arguments override the preset.
#' @param hidden_dim Dimension of each node's embedding.
#' @param weights `"edge"` (one learnable weight per edge, as in BINNs) or
#'   `"shared"` (fixed mean aggregation with learnable per-level transforms).
#' @param readout Which nodes feed the prediction head: `"roots"` (top-level
#'   pathways, as in BINNs) or `"pathways"` (all pathway nodes, akin to
#'   deep supervision / graph-level pooling).
#' @param activation Activation function: `"tanh"`, `"relu"`, `"gelu"`,
#'   `"sigmoid"` or `"identity"`.
#' @param dropout Dropout probability applied to hidden node embeddings.
#' @param gene_layer Whether to add a gene node combining all omics features
#'   of each gene before the first pathway layer (as in P-NET). `NULL`
#'   (default) uses a gene layer only for multi-omics inputs.
#' @return An object of class `bn_arch`.
#' @examples
#' bn_arch()                      # classic BINN
#' bn_arch("gnn")                 # relaxed graph neural network
#' bn_arch(hidden_dim = 4)        # BINN with 4-d pathway embeddings
#' @export
bn_arch <- function(preset = c("binn", "gnn"),
                    hidden_dim = NULL,
                    weights = NULL,
                    readout = NULL,
                    activation = NULL,
                    dropout = NULL,
                    gene_layer = NULL) {
  preset <- rlang::arg_match(preset)
  defaults <- switch(preset,
    binn = list(hidden_dim = 1L, weights = "edge", readout = "roots",
                activation = "tanh", dropout = 0.1),
    gnn = list(hidden_dim = 4L, weights = "shared", readout = "roots",
               activation = "tanh", dropout = 0.1)
  )
  out <- list(
    preset = preset,
    hidden_dim = as.integer(hidden_dim %||% defaults$hidden_dim),
    weights = rlang::arg_match0(weights %||% defaults$weights, c("edge", "shared")),
    readout = rlang::arg_match0(readout %||% defaults$readout, c("roots", "pathways")),
    activation = rlang::arg_match0(activation %||% defaults$activation,
                                   c("tanh", "relu", "gelu", "sigmoid", "identity")),
    dropout = dropout %||% defaults$dropout,
    gene_layer = gene_layer
  )
  if (length(out$hidden_dim) != 1 || is.na(out$hidden_dim) || out$hidden_dim < 1) {
    cli::cli_abort("{.arg hidden_dim} must be a positive integer.")
  }
  if (!is.numeric(out$dropout) || out$dropout < 0 || out$dropout >= 1) {
    cli::cli_abort("{.arg dropout} must be in [0, 1).")
  }
  structure(out, class = "bn_arch")
}

#' @export
print.bn_arch <- function(x, ...) {
  cat_line("<bn_arch> {x$preset} preset")
  cat_bullets(c(
    "node embedding dimension: {x$hidden_dim}",
    "edge weights: {x$weights}",
    "readout: {x$readout}",
    "activation: {x$activation}; dropout: {x$dropout}",
    "gene layer: {if (is.null(x$gene_layer)) 'auto' else x$gene_layer}"
  ))
  invisible(x)
}

#' Training options
#'
#' @param epochs Maximum number of passes through the training data.
#' @param batch_size Mini-batch size. `NULL` uses 64 for classification and
#'   regression, and the full training set for Cox models (whose partial
#'   likelihood depends on the risk set).
#' @param lr Learning rate for the AdamW optimiser.
#' @param weight_decay Decoupled weight decay (L2 regularisation).
#' @param l1 L1 penalty on edge weights (encourages sparse explanations).
#' @param validation Fraction of samples held out for early stopping
#'   (stratified by class for classification). Use `0` to disable.
#' @param patience Stop after this many epochs without improvement in
#'   validation loss; the best weights are restored.
#' @param scale Standardise omics features and numeric covariates using
#'   training-set means and standard deviations.
#' @param seed Optional integer seed for reproducible initialisation, data
#'   splitting and batching.
#' @param verbose Print progress.
#' @param device Torch device, e.g. `"cpu"` or `"cuda"`.
#' @return An object of class `bn_control`.
#' @examples
#' bn_control(epochs = 50, seed = 1)
#' @export
bn_control <- function(epochs = 200L,
                       batch_size = NULL,
                       lr = 5e-3,
                       weight_decay = 1e-3,
                       l1 = 0,
                       validation = 0.2,
                       patience = 20L,
                       scale = TRUE,
                       seed = NULL,
                       verbose = FALSE,
                       device = "cpu") {
  if (validation < 0 || validation >= 1) {
    cli::cli_abort("{.arg validation} must be in [0, 1).")
  }
  if (!is.numeric(epochs) || epochs < 1) cli::cli_abort("{.arg epochs} must be a positive integer.")
  if (!is.numeric(patience) || patience < 1) cli::cli_abort("{.arg patience} must be a positive integer.")
  if (!is.numeric(lr) || lr <= 0) cli::cli_abort("{.arg lr} must be positive.")
  if (!is.null(batch_size) && (!is.numeric(batch_size) || batch_size < 2)) {
    cli::cli_abort("{.arg batch_size} must be at least 2.")
  }
  structure(
    list(epochs = as.integer(epochs), batch_size = batch_size, lr = lr,
         weight_decay = weight_decay, l1 = l1, validation = validation,
         patience = as.integer(patience), scale = scale, seed = seed,
         verbose = verbose, device = device),
    class = "bn_control"
  )
}

#' @export
print.bn_control <- function(x, ...) {
  cat_line("<bn_control>")
  cat_bullets(c(
    "epochs: {x$epochs} (patience {x$patience}, validation {x$validation})",
    "AdamW lr {x$lr}, weight decay {x$weight_decay}, L1 {x$l1}",
    "batch size: {x$batch_size %||% 'auto'}; device: {x$device}"
  ))
  invisible(x)
}
