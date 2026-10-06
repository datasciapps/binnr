#' Knowledge graphs for biologically informed networks
#'
#' A `bn_graph` is a directed acyclic graph whose edges point from genes to
#' the pathways they belong to, and from child pathways to their parent
#' pathways. It is an [igraph][igraph::igraph-package] object with node
#' attributes `name`, `type` (`"gene"`, `"pathway"` or `"dummy"`), `level`
#' (longest path from a gene; genes are level 0) and `label`, so every igraph
#' function works on it, and `tidygraph::as_tbl_graph()` turns it into a tidy
#' graph. [bn_nodes()] and [bn_edges()] return the node and edge tables as
#' tibbles. Because the prior knowledge is kept as *data* rather than
#' hard-coded into an architecture, graphs can be filtered, joined, rewired or
#' padded before fitting a model with [bn_fit()].
#'
#' `bn_graph()` is the low-level constructor: supply an edge list and,
#' optionally, which nodes are genes. Most users will want
#' [pathway_graph()] or [reactome_graph()] instead.
#'
#' @param edges A data frame with character columns `from` and `to`, or an
#'   igraph object.
#' @param genes Character vector naming the gene (input-level) nodes. Defaults
#'   to all nodes with no incoming edges.
#' @param labels Optional named character vector of human-readable node
#'   labels (e.g. pathway names), named by node id.
#'
#' @return An object of class `bn_graph` (and `igraph`).
#' @seealso [pathway_graph()], [reactome_graph()], [bn_prune()], [bn_rewire()],
#'   [bn_pad()]
#' @examples
#' g <- bn_graph(data.frame(
#'   from = c("G1", "G2", "G2", "P1"),
#'   to   = c("P1", "P1", "P2", "P2")
#' ))
#' g
#' bn_nodes(g)
#' igraph::topo_sort(g)
#' @export
bn_graph <- function(edges, genes = NULL, labels = NULL) {
  if (igraph::is_igraph(edges)) {
    e <- igraph::as_data_frame(edges, "edges")[c("from", "to")]
    old <- igraph::as_data_frame(edges, "vertices")
    if (is.null(labels) && "label" %in% names(old)) {
      labels <- stats::setNames(old$label, old$name)
    }
    if (is.null(genes) && "type" %in% names(old)) genes <- old$name[old$type == "gene"]
    dummies <- if ("type" %in% names(old)) old$name[old$type == "dummy"] else character()
    edges <- e
  } else {
    dummies <- character()
  }
  if (!is.data.frame(edges) || !all(c("from", "to") %in% names(edges))) {
    cli::cli_abort("{.arg edges} must be a data frame with columns {.field from} and {.field to}.")
  }
  edges <- data.frame(from = as.character(edges$from), to = as.character(edges$to))
  if (anyNA(edges$from) || anyNA(edges$to)) {
    cli::cli_abort("{.arg edges} must not contain missing node names.")
  }
  edges <- unique(edges[edges$from != edges$to, ])
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  if (!igraph::is_dag(g)) {
    cyc <- igraph::V(g)$name[igraph::feedback_arc_set(g, algo = "approx_eades") |>
                               igraph::head_of(graph = g) |> as.integer()]
    cli::cli_abort(c(
      "The graph contains a cycle, so it is not a DAG.",
      i = "Nodes involved include {.val {utils::head(unique(cyc), 5)}}."
    ))
  }
  nodes <- igraph::V(g)$name
  if (is.null(genes)) {
    genes <- nodes[igraph::degree(g, mode = "in") == 0]
  } else {
    genes <- intersect(as.character(genes), nodes)
    if (any(igraph::degree(g, v = genes, mode = "in") > 0)) {
      cli::cli_abort("Gene nodes must not have incoming edges.")
    }
  }
  type <- ifelse(nodes %in% genes, "gene", "pathway")
  type[nodes %in% dummies] <- "dummy"
  igraph::V(g)$type <- type
  igraph::V(g)$level <- dag_levels(g)
  lab <- nodes
  if (!is.null(labels)) {
    hit <- nodes %in% names(labels)
    lab[hit] <- unname(labels[nodes[hit]])
  }
  igraph::V(g)$label <- lab
  # Stable node order: level, then type, then name
  ord <- order(igraph::V(g)$level, type, nodes)
  g <- igraph::permute(g, match(seq_along(nodes), ord))
  class(g) <- c("bn_graph", class(g))
  g
}

#' Longest-path level of each node in a DAG (sources have level 0)
#' @noRd
dag_levels <- function(g) {
  ord <- as.integer(igraph::topo_sort(g, mode = "out"))
  level <- integer(igraph::vcount(g))
  adj <- igraph::adjacent_vertices(g, ord, mode = "in")
  for (i in seq_along(ord)) {
    p <- as.integer(adj[[i]])
    if (length(p)) level[ord[i]] <- max(level[p]) + 1L
  }
  level
}

#' @rdname bn_graph
#' @param graph A `bn_graph`.
#' @export
bn_nodes <- function(graph) {
  check_graph(graph)
  tibble::as_tibble(igraph::as_data_frame(graph, "vertices"))[c("name", "type", "level", "label")]
}

#' @rdname bn_graph
#' @export
bn_edges <- function(graph) {
  check_graph(graph)
  tibble::as_tibble(igraph::as_data_frame(graph, "edges"))[c("from", "to")]
}

#' @export
print.bn_graph <- function(x, ...) {
  n <- bn_nodes(x)
  n_gene <- sum(n$type == "gene")
  n_path <- sum(n$type == "pathway")
  n_dummy <- sum(n$type == "dummy")
  n_edge <- igraph::ecount(x)
  n_root <- length(graph_roots(x))
  cat_line("<bn_graph> {n_gene} gene{?s}, {n_path} pathway{?s}, {n_edge} edge{?s}")
  cat_line("Depth {max(n$level)}; {n_root} root pathway{?s}",
           if (n_dummy > 0) "; {n_dummy} dummy (padding) node{?s}" else "")
  invisible(x)
}

#' @export
summary.bn_graph <- function(object, ...) {
  n <- bn_nodes(object)
  e <- bn_edges(object)
  lv <- n$level[match(e$to, n$name)] - n$level[match(e$from, n$name)]
  genes <- n$name[n$type == "gene"]
  tibble::tibble(
    n_genes = length(genes),
    n_pathways = sum(n$type == "pathway"),
    n_dummy = sum(n$type == "dummy"),
    n_edges = nrow(e),
    depth = max(n$level),
    n_roots = length(graph_roots(object)),
    skip_edges = sum(lv > 1L),
    density = sum(e$from %in% genes) / max(1, length(genes) * sum(n$type != "gene"))
  )
}

#' @rdname bn_graph
#' @param x A `bn_graph`.
#' @param ... Unused.
#' @details `as_tibble()` on a graph returns its edge table; `as.igraph()`
#'   drops the `bn_graph` class.
#' @exportS3Method tibble::as_tibble
as_tibble.bn_graph <- function(x, ...) {
  bn_edges(x)
}

#' @rdname bn_graph
#' @exportS3Method igraph::as.igraph
as.igraph.bn_graph <- function(x, ...) {
  class(x) <- setdiff(class(x), "bn_graph")
  x
}

graph_roots <- function(graph) {
  igraph::V(graph)$name[igraph::V(graph)$type != "gene" & igraph::degree(graph, mode = "out") == 0]
}

check_graph <- function(graph, arg = rlang::caller_arg(graph),
                        call = rlang::caller_env()) {
  if (!inherits(graph, "bn_graph")) {
    cli::cli_abort("{.arg {arg}} must be a {.cls bn_graph}, not {.obj_type_friendly {graph}}.",
                   call = call)
  }
  invisible(graph)
}

# Genes reachable to each non-gene node: a logical genes x pathways matrix
descendant_genes <- function(graph) {
  n <- bn_nodes(graph)
  genes <- n$name[n$type == "gene"]
  pws <- n$name[n$type != "gene"]
  if (!length(genes) || !length(pws)) {
    return(matrix(FALSE, length(genes), length(pws), dimnames = list(genes, pws)))
  }
  is.finite(igraph::distances(graph, v = genes, to = pws, mode = "out"))
}

#' Build a knowledge graph from gene sets and a pathway hierarchy
#'
#' Converts pathway membership lists (e.g. from Reactome, KEGG, GO or
#' `msigdbr::msigdbr()`) and an optional parent/child pathway hierarchy into a
#' [bn_graph].
#'
#' Many gene-set resources (including Reactome's GMT files and MSigDB) list
#' every gene annotated to a pathway *or any of its sub-pathways*. With
#' `redundant = "drop"` (the default, as in P-NET and BINN), a gene is linked
#' only to the most specific pathways containing it; information still reaches
#' ancestors through the hierarchy. With `redundant = "keep"`, genes are also
#' wired directly to every ancestor: a graph formulation handles these "skip"
#' edges natively, whereas a stratified BINN would need padding nodes.
#'
#' @param gene_sets Either a named list of character vectors (pathway ->
#'   genes) or a data frame in long format with one row per gene-pathway pair.
#' @param hierarchy Optional data frame with one row per parent/child pathway
#'   relation.
#' @param pathway_col,gene_col Column names in `gene_sets` (data-frame form).
#'   By default common names such as `pathway`, `gs_name`, `gene`,
#'   `gene_symbol` are detected.
#' @param parent_col,child_col Column names in `hierarchy`.
#' @param labels Optional named character vector (or two-column data frame of
#'   id and name) of pathway labels.
#' @param min_genes,max_genes Keep pathways annotated (at any level) to between
#'   `min_genes` and `max_genes` genes.
#' @param redundant `"drop"` or `"keep"` gene-to-ancestor edges; see Details.
#'
#' @return A [bn_graph].
#' @seealso [reactome_graph()], [bn_prune()]
#' @examples
#' sets <- list(P1 = c("A", "B"), P2 = c("B", "C"), P3 = c("A", "B", "C"))
#' hier <- data.frame(parent = c("P3", "P3"), child = c("P1", "P2"))
#' pathway_graph(sets, hier)
#'
#' # Built-in Reactome snapshot restricted to the genes of the TCGA example cohorts
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
#'                    labels = reactome_tcga$pathways)
#' g
#' @export
pathway_graph <- function(gene_sets,
                          hierarchy = NULL,
                          pathway_col = NULL,
                          gene_col = NULL,
                          parent_col = "parent",
                          child_col = "child",
                          labels = NULL,
                          min_genes = 1,
                          max_genes = Inf,
                          redundant = c("drop", "keep")) {
  redundant <- rlang::arg_match(redundant)
  memb <- as_membership(gene_sets, pathway_col, gene_col)
  if (!is.null(hierarchy)) {
    if (!all(c(parent_col, child_col) %in% names(hierarchy))) {
      cli::cli_abort("{.arg hierarchy} must have columns {.field {parent_col}} and {.field {child_col}}.")
    }
    hier <- unique(data.frame(from = as.character(hierarchy[[child_col]]),
                              to = as.character(hierarchy[[parent_col]])))
  } else {
    hier <- data.frame(from = character(), to = character())
  }
  # Graph with all-level membership, to count descendant genes per pathway
  full <- bn_graph(rbind(data.frame(from = memb$gene, to = memb$pathway), hier),
                   genes = unique(memb$gene))
  size <- colSums(descendant_genes(full))
  keep_pw <- names(size)[size >= min_genes & size <= max_genes]
  memb <- memb[memb$pathway %in% keep_pw, ]
  hier <- hier[hier$from %in% keep_pw & hier$to %in% keep_pw, ]

  if (redundant == "drop" && nrow(hier) > 0) {
    # A gene-pathway edge is redundant if the gene also reaches the pathway
    # through a child pathway, i.e. there is a path of length two.
    g2 <- bn_graph(rbind(data.frame(from = memb$gene, to = memb$pathway), hier),
                   genes = unique(memb$gene))
    A <- igraph::as_adjacency_matrix(g2, sparse = TRUE)
    two <- A %*% A
    via_child <- two[cbind(match(memb$gene, rownames(two)), match(memb$pathway, colnames(two)))] > 0
    memb <- memb[!via_child, ]
  } else if (redundant == "keep") {
    D <- descendant_genes(full)
    idx <- which(D[, keep_pw, drop = FALSE], arr.ind = TRUE)
    memb <- data.frame(pathway = keep_pw[idx[, 2]], gene = rownames(D)[idx[, 1]])
  }

  edges <- rbind(data.frame(from = memb$gene, to = memb$pathway), hier)
  bn_graph(edges, genes = unique(memb$gene), labels = as_labels(labels))
}

as_membership <- function(gene_sets, pathway_col = NULL, gene_col = NULL,
                          call = rlang::caller_env()) {
  if (is.list(gene_sets) && !is.data.frame(gene_sets)) {
    if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
      cli::cli_abort("A list of gene sets must be named by pathway.", call = call)
    }
    out <- data.frame(
      pathway = rep(names(gene_sets), lengths(gene_sets)),
      gene = as.character(unlist(gene_sets, use.names = FALSE))
    )
  } else if (is.data.frame(gene_sets)) {
    pathway_col <- pathway_col %||%
      pick_col(gene_sets, c("pathway", "gs_name", "gs_id", "term", "gene_set", "set"))
    gene_col <- gene_col %||%
      pick_col(gene_sets, c("gene", "gene_symbol", "symbol", "gene_id", "ensembl_gene", "feature"))
    if (is.null(pathway_col) || is.null(gene_col)) {
      cli::cli_abort(c(
        "Could not identify pathway and gene columns in {.arg gene_sets}.",
        i = "Specify {.arg pathway_col} and {.arg gene_col}."
      ), call = call)
    }
    out <- data.frame(
      pathway = as.character(gene_sets[[pathway_col]]),
      gene = as.character(gene_sets[[gene_col]])
    )
  } else {
    cli::cli_abort("{.arg gene_sets} must be a named list or a data frame.", call = call)
  }
  out <- out[!is.na(out$gene) & out$gene != "", ]
  unique(out)
}

pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) NULL else hit[[1]]
}

as_labels <- function(labels) {
  if (is.null(labels)) return(NULL)
  if (is.data.frame(labels)) {
    return(stats::setNames(as.character(labels[[2]]), as.character(labels[[1]])))
  }
  labels
}

#' Restrict a knowledge graph to observed genes
#'
#' Removes gene nodes not in `genes`, then removes pathways that no longer
#' have any gene among their descendants, and (optionally) pathways with fewer
#' than `min_genes` descendant genes or above a maximum level. [bn_fit()]
#' calls this automatically with the genes present in the data.
#'
#' @param graph A [bn_graph].
#' @param genes Character vector of genes to keep. Defaults to all genes.
#' @param min_genes Minimum number of descendant genes for a pathway to be
#'   kept.
#' @param max_level Optional maximum pathway level; higher pathways are
#'   removed (so that nodes at `max_level` become roots). Genes whose only
#'   annotations lie above `max_level` are dropped too. Feedforward BINNs such
#'   as P-NET fix this depth (e.g. 5); binnr uses the full hierarchy by default.
#' @return A [bn_graph].
#' @examplesIf binnr::tcga_cached("BRCA")
#' brca <- tcga_cohort("BRCA")
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
#' bn_prune(g, genes = colnames(brca$omics$rna)[1:100], min_genes = 3)
#' bn_prune(g, max_level = 5)
#' @export
bn_prune <- function(graph, genes = NULL, min_genes = 1, max_level = Inf) {
  check_graph(graph)
  n <- bn_nodes(graph)
  all_genes <- n$name[n$type == "gene"]
  genes <- if (is.null(genes)) all_genes else intersect(all_genes, genes)
  if (length(genes) == 0) {
    cli::cli_abort("None of the supplied {.arg genes} are in the graph.")
  }
  g <- igraph::delete_vertices(graph, setdiff(all_genes, genes))
  g <- igraph::delete_vertices(g, n$name[n$type != "gene" & n$level > max_level])
  class(g) <- c("bn_graph", setdiff(class(g), "bn_graph"))
  size <- colSums(descendant_genes(g))
  g <- igraph::delete_vertices(g, names(size)[size < min_genes])
  # drop genes left without any pathway
  g <- igraph::delete_vertices(g, igraph::V(g)$name[igraph::degree(g, mode = "out") == 0 &
                                                     igraph::V(g)$type == "gene"])
  bn_graph(g)
}

#' Download a Reactome knowledge graph
#'
#' Downloads the current Reactome gene sets and pathway hierarchy
#' (`ReactomePathwaysRelation.txt`) for a species, caches them, and returns a
#' [bn_graph] via [pathway_graph()]. With `identifiers = "symbol"` the gene
#' sets come from `ReactomePathways.gmt` (HGNC symbols, all pathway levels);
#' with `"uniprot"` they come from `UniProt2Reactome.txt` (UniProt accessions
#' annotated to their lowest-level pathways), which is what proteomics data
#' such as [hartman_aki()] use.
#'
#' @param species Species prefix used by Reactome stable identifiers, e.g.
#'   `"HSA"` (human) or `"MMU"` (mouse). The symbol GMT file is human only;
#'   the UniProt mapping covers all Reactome species.
#' @param identifiers `"symbol"` (gene symbols) or `"uniprot"` (UniProt
#'   accessions) for the input nodes.
#' @param cache_dir Directory where downloaded files are cached (see
#'   [binnr_cache_dir()]).
#' @param refresh Re-download even if cached files exist.
#' @param ... Passed to [pathway_graph()], e.g. `min_genes`, `max_genes`,
#'   `redundant`.
#' @return A [bn_graph] with pathway names as node labels.
#' @examplesIf interactive()
#' g <- reactome_graph(min_genes = 5)
#' g
#' reactome_graph(identifiers = "uniprot")
#' @export
reactome_graph <- function(species = "HSA", identifiers = c("symbol", "uniprot"),
                           cache_dir = binnr_cache_dir(), refresh = FALSE, ...) {
  identifiers <- rlang::arg_match(identifiers)
  base <- "https://reactome.org/download/current/"
  prefix <- paste0("R-", species, "-")
  rel_path <- cached_download(paste0(base, "ReactomePathwaysRelation.txt"),
                              "reactome/ReactomePathwaysRelation.txt", cache_dir, refresh)
  names_path <- cached_download(paste0(base, "ReactomePathways.txt"),
                                "reactome/ReactomePathways.txt", cache_dir, refresh)
  if (identifiers == "symbol") {
    gmt_path <- cached_download(paste0(base, "ReactomePathways.gmt.zip"),
                                "reactome/ReactomePathways.gmt.zip", cache_dir, refresh)
    gmt <- strsplit(readLines(unz(gmt_path, "ReactomePathways.gmt")), "\t", fixed = TRUE)
    gene_sets <- data.frame(
      pathway = rep(vapply(gmt, `[`, "", 2L), lengths(gmt) - 2L),
      gene = unlist(lapply(gmt, `[`, -(1:2)), use.names = FALSE)
    )
  } else {
    up_path <- cached_download(paste0(base, "UniProt2Reactome.txt"),
                               "reactome/UniProt2Reactome.txt", cache_dir, refresh)
    up <- utils::read.delim(up_path, header = FALSE, quote = "",
                            col.names = c("gene", "pathway", "url", "name", "evidence", "species"))
    gene_sets <- unique(up[startsWith(up$pathway, prefix), c("pathway", "gene")])
  }
  gene_sets <- gene_sets[startsWith(gene_sets$pathway, prefix), ]
  rel <- utils::read.delim(rel_path, header = FALSE, col.names = c("parent", "child"))
  rel <- rel[startsWith(rel$parent, prefix) & startsWith(rel$child, prefix), ]
  nm <- utils::read.delim(names_path, header = FALSE, quote = "",
                          col.names = c("pathway", "name", "species"))
  nm <- nm[startsWith(nm$pathway, prefix), c("pathway", "name")]
  pathway_graph(gene_sets, rel, labels = nm, ...)
}

cat_line <- function(..., .envir = parent.frame()) {
  cat(cli::format_inline(..., .envir = .envir), "\n", sep = "")
}

cat_bullets <- function(x, .envir = parent.frame()) {
  for (b in x) cat("* ", cli::format_inline(b, .envir = .envir), "\n", sep = "")
}
