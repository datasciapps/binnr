#' Randomised null knowledge graphs
#'
#' Creates a randomised version of a [bn_graph] with the same sparsity, to
#' test whether a model benefits from the *specific* prior knowledge or merely
#' from sparsity (cf. "sparsity is all you need"). Only gene-to-pathway edges
#' are randomised; the pathway hierarchy is kept. Null models differ in what
#' they preserve:
#'
#' * `"degree"` (default): degree-preserving rewiring (Maslov & Sneppen,
#'   2002). Pairs of gene-pathway edges `(g1, p1), (g2, p2)` are repeatedly
#'   swapped to `(g1, p2), (g2, p1)`, so every gene keeps its number of
#'   pathways and every pathway keeps its number of genes.
#' * `"pathway_size"`: each pathway keeps its number of genes, drawn uniformly
#'   from all genes; gene degrees are not preserved.
#' * `"bernoulli"`: an Erdos-Renyi graph with the same number of gene-pathway
#'   edges, placed uniformly at random.
#' * `"genes"`: permutes gene labels, so the topology is identical and only
#'   the assignment of genes to positions in the graph is random.
#'
#' What such a comparison can show is limited. A gene's influence on the
#' output is bounded by the paths it has to the readout, so a null model
#' that changes how many paths the *informative* genes keep is testing that,
#' not the biology. `"degree"` and `"genes"` preserve each gene's number of
#' direct pathway memberships but not its number of paths to the roots;
#' `"bernoulli"` and `"pathway_size"` preserve neither. Comparing several null
#' models, and checking that the model can exploit a *known* planted signal at
#' all, guards against over-interpreting a single "real vs random" result.
#' Genes that end up with no pathway are dropped from `"bernoulli"` and
#' `"pathway_size"` graphs.
#'
#' @param graph A [bn_graph].
#' @param method One of `"degree"`, `"pathway_size"`, `"bernoulli"`, `"genes"`.
#' @param n_swaps For `"degree"`, number of attempted swaps per edge.
#' @param seed Optional random seed (the global RNG state is restored).
#' @return A [bn_graph] with attribute `"null_model"`.
#' @examples
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
#' g_null <- bn_rewire(g, seed = 1)
#' # Same degree sequence, different edges
#' identical(sort(table(bn_edges(g)$from)), sort(table(bn_edges(g_null)$from)))
#' summary(bn_rewire(g, "bernoulli", seed = 1))
#' @export
bn_rewire <- function(graph, method = c("degree", "pathway_size", "bernoulli", "genes"),
                      n_swaps = 10, seed = NULL) {
  check_graph(graph)
  method <- rlang::arg_match(method)
  if (!is.null(seed)) withr::local_seed(seed)
  n <- bn_nodes(graph)
  e <- bn_edges(graph)
  genes <- n$name[n$type == "gene"]
  is_ge <- e$from %in% genes
  ge <- e[is_ge, ]
  pw <- unique(ge$to)
  if (method == "genes") {
    perm <- stats::setNames(sample(genes), genes)
    ge$from <- unname(perm[ge$from])
  } else if (method == "pathway_size") {
    size <- table(ge$to)[pw]
    ge <- tibble::tibble(
      from = unlist(lapply(size, function(k) sample(genes, k)), use.names = FALSE),
      to = rep(pw, times = size)
    )
  } else if (method == "bernoulli") {
    # G(n, m) bipartite graph with the same number of gene-pathway edges
    b <- igraph::sample_bipartite_gnm(length(genes), length(pw),
                                      m = nrow(ge), directed = TRUE, mode = "out")
    be <- igraph::as_edgelist(b, names = FALSE)
    ge <- tibble::tibble(from = genes[be[, 1]], to = pw[be[, 2] - length(genes)])
  } else {
    # Degree-preserving rewiring of the gene-pathway subgraph: igraph swaps
    # (g1, p1), (g2, p2) -> (g1, p2), (g2, p1), keeping in- and out-degrees
    sub <- igraph::graph_from_data_frame(ge, directed = TRUE)
    sub <- igraph::rewire(sub, igraph::keeping_degseq(niter = ceiling(n_swaps * nrow(ge))))
    ge <- tibble::as_tibble(igraph::as_data_frame(sub, "edges"))
  }
  out <- bn_graph(rbind(ge, e[!is_ge, ]), genes = intersect(genes, ge$from),
                  labels = stats::setNames(n$label, n$name))
  out <- mark_dummies(out, n$name[n$type == "dummy"])
  attr(out, "null_model") <- method
  out
}

#' Pad a knowledge graph into strata (as in feedforward BINNs)
#'
#' Feedforward BINN implementations update the network stratum by stratum
#' (a *stratified* schedule), so every edge must connect adjacent strata.
#' "Skip" edges (e.g. a gene annotated directly to a high-level pathway) and
#' roots at different depths are handled by inserting copy ("dummy") nodes,
#' as in P-NET and BINN. `bn_pad()` performs this transformation so that the
#' cost of artificial stratification can be inspected, or so that a
#' stratified BINN can be reproduced exactly (fit the padded graph with
#' `bn_arch("binn")`). Models in binnr do *not* need it: the *topological*
#' schedule updates each node once, after its children, and handles skip
#' edges natively. Note that a padded copy node with its own weight and bias
#' is not equivalent to a skip connection, so the two graphs define different
#' function classes.
#'
#' @param graph A [bn_graph].
#' @param roots If `TRUE`, also pad root pathways up to the maximum depth, so
#'   that all roots lie in the final layer.
#' @return A [bn_graph] whose extra nodes have `type == "dummy"`.
#' @examples
#' g <- bn_graph(data.frame(from = c("G1", "P1", "G2"), to = c("P1", "P2", "P2")))
#' bn_pad(g)
#' bn_edges(bn_pad(g))
#' @export
bn_pad <- function(graph, roots = TRUE) {
  check_graph(graph)
  n <- bn_nodes(graph)
  e <- bn_edges(graph)
  lvl <- stats::setNames(n$level, n$name)
  gap <- lvl[e$to] - lvl[e$from]
  new_edges <- list(e[gap == 1L, ])
  dummy <- character()
  long <- e[gap > 1L, ]
  copy_name <- function(u, l) paste0(u, "@", l)
  for (i in seq_len(nrow(long))) {
    u <- long$from[i]
    v <- long$to[i]
    lu <- lvl[[u]]
    lv <- lvl[[v]]
    chain <- c(u, copy_name(u, seq(lu + 1L, lv - 1L)), v)
    dummy <- c(dummy, chain[-c(1L, length(chain))])
    new_edges[[length(new_edges) + 1L]] <- tibble::tibble(
      from = chain[-length(chain)], to = chain[-1L]
    )
  }
  if (roots) {
    L <- max(lvl)
    for (r in graph_roots(graph)) {
      if (lvl[[r]] < L) {
        chain <- c(r, copy_name(r, seq(lvl[[r]] + 1L, L)))
        dummy <- c(dummy, chain[-1L])
        new_edges[[length(new_edges) + 1L]] <- tibble::tibble(
          from = chain[-length(chain)], to = chain[-1L]
        )
      }
    }
  }
  edges <- do.call(rbind, new_edges)
  labels <- stats::setNames(n$label, n$name)
  out <- bn_graph(edges, genes = n$name[n$type == "gene"], labels = labels)
  d <- unique(dummy)
  mark_dummies(out, d, labels = paste0(labels[sub("@[0-9]+$", "", d)], " (copy)"))
}

# Set node type "dummy" (and optionally labels) on padding nodes
mark_dummies <- function(graph, dummies, labels = NULL) {
  idx <- match(intersect(dummies, igraph::V(graph)$name), igraph::V(graph)$name)
  if (!length(idx)) return(graph)
  igraph::V(graph)$type[idx] <- "dummy"
  if (!is.null(labels)) {
    igraph::V(graph)$label[idx] <- labels[match(igraph::V(graph)$name[idx], dummies)]
  }
  graph
}

#' Simulate an incomplete knowledge graph
#'
#' Randomly removes a proportion of gene-to-pathway annotations, mimicking
#' missing or not-yet-curated prior knowledge. Genes left without any pathway
#' and pathways left without any gene are removed. Use it to study how
#' robust a model is to incomplete knowledge graphs, e.g. by comparing fits on
#' `bn_drop_edges(g, 0.5)` with fits on `g`.
#'
#' @param graph A [bn_graph].
#' @param prop Proportion of gene-pathway edges to remove, in `[0, 1)`.
#' @param seed Optional random seed (the global RNG state is restored).
#' @return A [bn_graph].
#' @examples
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
#' summary(g)
#' summary(bn_drop_edges(g, prop = 0.5, seed = 1))
#' @export
bn_drop_edges <- function(graph, prop = 0.5, seed = NULL) {
  check_graph(graph)
  if (!is.numeric(prop) || length(prop) != 1 || prop < 0 || prop >= 1) {
    cli::cli_abort("{.arg prop} must be a single number in [0, 1).")
  }
  if (!is.null(seed)) withr::local_seed(seed)
  n <- bn_nodes(graph)
  e <- bn_edges(graph)
  genes <- n$name[n$type == "gene"]
  ge <- which(e$from %in% genes)
  drop <- ge[sample.int(length(ge), round(prop * length(ge)))]
  kept <- if (length(drop)) e[-drop, ] else e
  out <- bn_graph(kept, genes = intersect(genes, kept$from),
                  labels = stats::setNames(n$label, n$name))
  out <- bn_prune(out)
  attr(out, "dropped") <- prop
  out
}
