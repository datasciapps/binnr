# The central claim: with scalar nodes, edge-specific weights and root readout,
# DAG message passing is numerically identical to a masked feedforward BINN.

masked_dense_forward <- function(fit, x) {
  p <- fit$params
  net <- fit$network
  nodes <- net$nodes
  genes <- sub("\\|x$", "", nodes$name[nodes$type == "input"])
  h <- scale(x[, genes], center = fit$scaler$center, scale = fit$scaler$scale)
  names_prev <- nodes$name[nodes$level == 0]
  for (k in seq_len(max(nodes$level))) {
    names_k <- nodes$name[nodes$level == k]
    # Dense weight matrix W (n_k x n_prev) with mask M from the graph
    W <- matrix(0, length(names_k), length(names_prev), dimnames = list(names_k, names_prev))
    e <- which(nodes$level[net$dst] == k)
    W[cbind(nodes$name[net$dst[e]], nodes$name[net$src[e]])] <- p$edge_weight[e]
    b <- p$bias[match(names_k, nodes$name) - net$n_input, 1]
    h <- tanh(h %*% t(W) + rep(b, each = nrow(h)))
    names_prev <- names_k
  }
  drop(h %*% t(p$head.weight) + p$head.bias[1])
}

test_that("BINN preset equals a masked dense feedforward network (layered graph)", {
  skip_if_no_torch()
  d <- toy_data()
  fit <- quiet_fit(d$x, d$y, graph = toy_layered(),
                   arch = bn_arch("binn", dropout = 0),
                   control = bn_control(epochs = 3, validation = 0, seed = 1))
  mp <- binnr:::predict_raw(fit, d$x)[, 1]
  ref <- masked_dense_forward(fit, d$x)
  expect_equal(mp, ref, tolerance = 1e-5)
})

test_that("message passing on a ragged DAG matches a manual topological computation", {
  skip_if_no_torch()
  d <- toy_data()
  fit <- quiet_fit(d$x, d$y, graph = toy_ragged(),
                   arch = bn_arch("binn", dropout = 0),
                   control = bn_control(epochs = 3, validation = 0, seed = 2))
  p <- fit$params
  net <- fit$network
  nodes <- net$nodes
  xs <- scale(d$x, center = fit$scaler$center, scale = fit$scaler$scale)
  H <- matrix(0, nrow(xs), nrow(nodes), dimnames = list(NULL, nodes$name))
  H[, seq_len(net$n_input)] <- xs
  for (v in which(nodes$level > 0)[order(nodes$level[nodes$level > 0])]) {
    e <- which(net$dst == v)
    z <- H[, net$src[e], drop = FALSE] %*% p$edge_weight[e] + p$bias[v - net$n_input, 1]
    H[, v] <- tanh(z)
  }
  ref <- drop(H[, sort(net$roots), drop = FALSE] %*% t(p$head.weight) + p$head.bias[1])
  expect_equal(binnr:::predict_raw(fit, d$x)[, 1], ref, tolerance = 1e-5)
  # skip edge G3 -> P3 is used directly: no dummy nodes in the network
  expect_false(any(nodes$type == "dummy"))
})

test_that("shared weights implement mean aggregation", {
  skip_if_no_torch()
  d <- toy_data()
  fit <- quiet_fit(d$x, d$y, graph = toy_layered(),
                   arch = bn_arch("binn", weights = "shared", dropout = 0),
                   control = bn_control(epochs = 2, validation = 0, seed = 3))
  expect_null(fit$params$edge_weight)
  expect_equal(dim(fit$params$level_weight), c(2L, 1L, 1L))
  p <- fit$params
  xs <- scale(d$x, center = fit$scaler$center, scale = fit$scaler$scale)
  w1 <- p$level_weight[1, 1, 1]; w2 <- p$level_weight[2, 1, 1]
  P1 <- tanh(w1 * rowMeans(xs[, c("G1", "G2")]) + p$bias[1, 1])
  P2 <- tanh(w1 * rowMeans(xs[, c("G2", "G3", "G4")]) + p$bias[2, 1])
  P3 <- tanh(w2 * (P1 + P2) / 2 + p$bias[3, 1])
  expect_equal(binnr:::predict_raw(fit, d$x)[, 1],
               P3 * p$head.weight[1, 1] + p$head.bias[1], tolerance = 1e-5)
})
