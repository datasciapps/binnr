# Compile a knowledge graph plus observed features into the integer message
# passing schedule used by the torch module.
compile_network <- function(graph, features, gene_layer) {
  modalities <- names(features)
  genes_obs <- unique(unlist(features, use.names = FALSE))
  g <- bn_prune(graph, genes = genes_obs)
  gn <- bn_nodes(g)
  ge <- bn_edges(g)
  graph_genes <- gn$name[gn$type == "gene"]

  keep <- lapply(features, function(f) which(f %in% graph_genes))
  inputs <- tibble::tibble(
    modality = rep(modalities, lengths(keep)),
    gene = unlist(Map(function(f, k) f[k], features, keep), use.names = FALSE)
  )
  inputs$name <- paste(inputs$gene, inputs$modality, sep = "|")

  is_gene_edge <- ge$from %in% graph_genes
  if (gene_layer) {
    edges <- rbind(
      tibble::tibble(from = inputs$name, to = inputs$gene),
      ge
    )
    types <- c(stats::setNames(rep("input", nrow(inputs)), inputs$name),
               stats::setNames(gn$type, gn$name))
  } else {
    gpe <- ge[is_gene_edge, ]
    m <- merge(inputs[, c("name", "gene")], gpe, by.x = "gene", by.y = "from")
    edges <- rbind(
      tibble::tibble(from = m$name, to = m$to),
      ge[!is_gene_edge, ]
    )
    gn2 <- gn[gn$type != "gene", ]
    types <- c(stats::setNames(rep("input", nrow(inputs)), inputs$name),
               stats::setNames(gn2$type, gn2$name))
  }
  node_names <- unique(c(inputs$name, edges$to, edges$from))
  ig <- igraph::graph_from_data_frame(edges, directed = TRUE,
                                      vertices = data.frame(name = node_names))
  lvl <- dag_levels(ig)
  lvl[match(inputs$name, node_names)] <- 0L
  # Stable order: level, then inputs in data order
  inp_pos <- match(node_names, inputs$name)
  ord <- order(lvl, ifelse(is.na(inp_pos), Inf, inp_pos), node_names)
  node_names <- node_names[ord]
  lvl <- lvl[ord]
  labels <- stats::setNames(gn$label, gn$name)
  nodes <- tibble::tibble(
    name = node_names,
    type = unname(types[node_names]),
    level = lvl,
    label = ifelse(node_names %in% names(labels), labels[node_names], node_names)
  )
  nodes$modality <- inputs$modality[match(nodes$name, inputs$name)]
  nodes$gene <- inputs$gene[match(nodes$name, inputs$name)]

  src <- match(edges$from, node_names)
  dst <- match(edges$to, node_names)
  indeg <- tabulate(dst, nbins = length(node_names))
  n_levels <- max(lvl)
  # Position of each node within its own level (levels are contiguous blocks)
  pos_in_level <- stats::ave(seq_along(lvl), lvl, FUN = seq_along)
  schedule <- lapply(seq_len(n_levels), function(k) {
    nodes_k <- which(lvl == k)
    e_k <- which(lvl[dst] == k)
    # Group incoming edges by the level of their source, so messages can be
    # gathered from each earlier level without concatenating node states.
    e_k <- e_k[order(lvl[src[e_k]], e_k)]
    src_lvl <- lvl[src[e_k]]
    groups <- lapply(sort(unique(src_lvl)), function(j) {
      list(level = j, src_local = pos_in_level[src[e_k[src_lvl == j]]])
    })
    list(
      nodes = nodes_k,
      edge_id = e_k,
      src = src[e_k],
      groups = groups,
      dst_local = match(dst[e_k], nodes_k),
      norm = 1 / indeg[dst[e_k]]
    )
  })
  n_input <- nrow(inputs)
  is_out <- nodes$type != "input" & !seq_along(node_names) %in% src
  hidden_pw <- which(nodes$type %in% c("pathway", "dummy"))
  structure(list(
    nodes = nodes,
    edges = tibble::tibble(from = edges$from, to = edges$to),
    src = src, dst = dst, indeg = indeg,
    schedule = schedule,
    n_input = n_input,
    input_modality = match(inputs$modality, modalities),
    modalities = modalities,
    keep = keep,
    roots = which(is_out),
    pathways = hidden_pw,
    gene_layer = gene_layer
  ), class = "bn_network")
}

activation_fn <- function(name) {
  switch(name,
    tanh = torch::torch_tanh,
    relu = torch::nnf_relu,
    gelu = torch::nnf_gelu,
    sigmoid = torch::torch_sigmoid,
    identity = function(x) x
  )
}

# The message-passing module ---------------------------------------------------
bn_module <- torch::nn_module(
  "bn_module",
  initialize = function(net, arch, n_out, n_cov, out_bias = TRUE) {
    d <- arch$hidden_dim
    self$d <- d
    self$n_input <- net$n_input
    self$n_nodes <- nrow(net$nodes)
    self$use_edge <- arch$weights == "edge"
    self$use_transform <- d > 1 || !self$use_edge
    self$act <- activation_fn(arch$activation)
    self$drop <- torch::nn_dropout(arch$dropout)
    long <- function(x) torch::torch_tensor(as.integer(x), dtype = torch::torch_long())
    flt <- function(x) torch::torch_tensor(as.numeric(x), dtype = torch::torch_float())

    n_levels <- length(net$schedule)
    self$n_levels <- n_levels
    # Buffers must be registered individually to move with $to(device)
    self$groups <- vector("list", n_levels)
    for (k in seq_len(n_levels)) {
      s <- net$schedule[[k]]
      self$groups[[k]] <- vapply(s$groups, function(gr) gr$level, integer(1))
      for (gi in seq_along(s$groups)) {
        self$register_buffer(paste0("src_", k, "_", gi), long(s$groups[[gi]]$src_local))
      }
      self$register_buffer(paste0("dst_", k), long(s$dst_local))
      self$register_buffer(paste0("eid_", k), long(s$edge_id))
      self$register_buffer(paste0("norm_", k), flt(s$norm))
      self$register_buffer(paste0("bidx_", k), long(s$nodes - net$n_input))
    }
    self$n_per_level <- vapply(net$schedule, function(s) length(s$nodes), integer(1))
    readout <- sort(if (arch$readout == "roots") net$roots else union(net$roots, net$pathways))
    self$register_buffer("readout", long(readout))
    # readout nodes grouped by level, for gathering without a full concatenation
    lvl <- net$nodes$level
    pos <- stats::ave(seq_along(lvl), lvl, FUN = seq_along)
    self$readout_levels <- sort(unique(lvl[readout]))
    for (j in self$readout_levels) {
      self$register_buffer(paste0("read_", j), long(pos[readout[lvl[readout] == j]]))
    }
    self$register_buffer("input_mod", long(net$input_modality))

    n_edges <- length(net$src)
    if (self$use_edge) {
      init_sd <- 1 / sqrt(pmax(net$indeg[net$dst], 1))
      self$edge_weight <- torch::nn_parameter(flt(stats::rnorm(n_edges, 0, init_sd)))
    }
    self$bias <- torch::nn_parameter(torch::torch_zeros(self$n_nodes - net$n_input, d))
    if (self$use_transform) {
      self$level_weight <- torch::nn_parameter(
        torch::torch_randn(n_levels, d, d) / sqrt(d)
      )
    }
    if (d > 1) {
      self$embed <- torch::nn_parameter(torch::torch_randn(length(net$modalities), d))
    }
    n_read <- length(readout) * d
    self$head <- torch::nn_linear(n_read + n_cov, n_out, bias = out_bias)
    self$n_cov <- n_cov
  },

  propagate = function(x) {
    # x: batch x n_input
    h0 <- x$unsqueeze(3)
    if (self$d > 1) {
      h0 <- h0 * self$embed$index_select(1, self$input_mod)$unsqueeze(1)
    }
    levels <- list(h0)
    for (k in seq_len(self$n_levels)) {
      src_levels <- self$groups[[k]]
      parts <- lapply(seq_along(src_levels), function(gi) {
        levels[[src_levels[[gi]] + 1L]]$index_select(2, self[[paste0("src_", k, "_", gi)]])
      })
      m <- if (length(parts) == 1L) parts[[1]] else torch::torch_cat(parts, dim = 2)
      if (self$use_transform) {
        m <- m$matmul(self$level_weight[k, , ])
      }
      w <- if (self$use_edge) {
        self$edge_weight$index_select(1, self[[paste0("eid_", k)]])
      } else {
        self[[paste0("norm_", k)]]
      }
      m <- m * w$view(c(1, -1, 1))
      agg <- torch::torch_zeros(x$size(1), self$n_per_level[[k]], self$d,
                                device = x$device, dtype = m$dtype)
      agg <- agg$index_add(2, self[[paste0("dst_", k)]], m)
      b <- self$bias$index_select(1, self[[paste0("bidx_", k)]])
      h <- self$drop(self$act(agg + b$unsqueeze(1)))
      levels[[k + 1L]] <- h
    }
    levels
  },

  readout_head = function(levels, cov = NULL) {
    parts <- lapply(self$readout_levels, function(j) {
      levels[[j + 1L]]$index_select(2, self[[paste0("read_", j)]])
    })
    z <- if (length(parts) == 1L) parts[[1]] else torch::torch_cat(parts, dim = 2)
    z <- z$flatten(start_dim = 2)
    if (!is.null(cov)) z <- torch::torch_cat(list(z, cov), dim = 2)
    self$head(z)
  },

  forward = function(x, cov = NULL) {
    self$readout_head(self$propagate(x), cov)
  }
)
