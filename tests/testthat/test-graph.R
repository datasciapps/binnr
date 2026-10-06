test_that("bn_graph computes longest-path levels and types", {
  g <- toy_ragged()
  expect_s3_class(g, "bn_graph")
  n <- bn_nodes(g)
  lvl <- stats::setNames(n$level, n$name)
  expect_equal(unname(lvl[c("G1", "P1", "P2", "P3", "P4")]), c(0L, 1L, 1L, 2L, 2L))
  expect_setequal(n$name[n$type == "gene"], paste0("G", 1:4))
  expect_setequal(binnr:::graph_roots(g), c("P3", "P4"))
  expect_equal(summary(g)$skip_edges, 2L) # G3 -> P3 and G4 -> P4
  expect_output(print(g), "4 genes")
})

test_that("cycles and malformed edges are rejected", {
  expect_error(bn_graph(data.frame(from = c("A", "B"), to = c("B", "A"))), "cycle")
  expect_error(bn_graph(data.frame(a = 1)), "from")
  expect_error(bn_graph(data.frame(from = c("G", "P"), to = c("P", "G")), genes = "G"),
               "incoming|cycle")
})

test_that("pathway_graph drops redundant ancestor edges by default", {
  sets <- list(P1 = c("A", "B"), P2 = c("B", "C"), P3 = c("A", "B", "C", "D"))
  hier <- data.frame(parent = c("P3", "P3"), child = c("P1", "P2"))
  g <- pathway_graph(sets, hier)
  e <- bn_edges(g)
  # A, B, C reach P3 via children; only D links to P3 directly
  expect_setequal(e$from[e$to == "P3"], c("P1", "P2", "D"))
  gk <- pathway_graph(sets, hier, redundant = "keep")
  expect_setequal(bn_edges(gk)$from[bn_edges(gk)$to == "P3"], c("P1", "P2", "A", "B", "C", "D"))
})

test_that("pathway_graph accepts long data frames (e.g. msigdbr)", {
  df <- data.frame(gs_name = c("S1", "S1", "S2"), gene_symbol = c("A", "B", "B"))
  g <- pathway_graph(df, min_genes = 2)
  expect_equal(nrow(bn_nodes(g)[bn_nodes(g)$type == "pathway", ]), 1)
  expect_error(pathway_graph(data.frame(x = 1, y = 2)), "identify")
})

test_that("pathway labels are attached", {
  g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
                     labels = reactome_tcga$pathways)
  n <- bn_nodes(g)
  expect_true(any(n$label == "Signal Transduction"))
})

test_that("bn_prune removes unobserved genes and empty pathways", {
  g <- toy_ragged()
  p <- bn_prune(g, genes = c("G1", "G2"))
  expect_setequal(bn_nodes(p)$name, c("G1", "G2", "P1", "P2", "P3", "P4"))
  p2 <- bn_prune(g, genes = c("G1", "G2", "G3", "G4"), min_genes = 3)
  expect_setequal(bn_nodes(p2)$name[bn_nodes(p2)$type == "pathway"], "P3")
  expect_error(bn_prune(g, genes = "nope"), "None")
})

test_that("bn_rewire preserves degrees and sparsity", {
  g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
  r <- bn_rewire(g, seed = 42)
  e0 <- bn_edges(g); e1 <- bn_edges(r)
  genes <- bn_nodes(g)$name[bn_nodes(g)$type == "gene"]
  ge0 <- e0[e0$from %in% genes, ]; ge1 <- e1[e1$from %in% genes, ]
  expect_equal(nrow(e0), nrow(e1))
  expect_equal(as.vector(table(ge0$from)[genes]), as.vector(table(ge1$from)[genes]))
  pw <- unique(ge0$to)
  expect_equal(as.vector(table(ge0$to)[pw]), as.vector(table(ge1$to)[pw]))
  overlap <- mean(paste(ge1$from, ge1$to) %in% paste(ge0$from, ge0$to))
  expect_lt(overlap, 0.5)
  # hierarchy is untouched
  expect_equal(e0[!e0$from %in% genes, ], e1[!e1$from %in% genes, ])
  # reproducible
  expect_equal(bn_edges(bn_rewire(g, seed = 42)), e1)
  gp <- bn_rewire(g, method = "genes", seed = 1)
  expect_equal(nrow(bn_edges(gp)), nrow(e0))
})

test_that("bn_pad produces a strictly layered graph with dummy nodes", {
  g <- toy_ragged()
  p <- bn_pad(g)
  s <- summary(p)
  expect_equal(s$skip_edges, 0L)
  expect_true(s$n_dummy >= 2)
  n <- bn_nodes(p)
  # all roots at the maximum level
  roots <- binnr:::graph_roots(p)
  expect_true(all(n$level[match(roots, n$name)] == max(n$level)))
  # Reactome example: padding inflates the graph
  gr <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
  expect_gt(summary(bn_pad(gr))$n_dummy, 1000)
})

test_that("bn_rewire keeps dummy nodes of padded graphs", {
  p <- bn_pad(toy_ragged())
  r <- bn_rewire(p, seed = 1)
  expect_equal(sum(bn_nodes(r)$type == "dummy"), sum(bn_nodes(p)$type == "dummy"))
})

test_that("bn_drop_edges removes the requested share of annotations", {
  g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
  genes <- bn_nodes(g)$name[bn_nodes(g)$type == "gene"]
  n0 <- sum(bn_edges(g)$from %in% genes)
  d <- bn_drop_edges(g, 0.5, seed = 1)
  n1 <- sum(bn_edges(d)$from %in% bn_nodes(d)$name[bn_nodes(d)$type == "gene"])
  expect_equal(n1, round(n0 * 0.5), tolerance = 1)
  expect_equal(bn_edges(bn_drop_edges(g, 0.5, seed = 1)), bn_edges(d))
  expect_equal(nrow(bn_edges(bn_drop_edges(g, 0))), nrow(bn_edges(g)))
  expect_error(bn_drop_edges(g, 1), "prop")
})

test_that("bernoulli and pathway_size null models preserve edge counts", {
  g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
  genes <- bn_nodes(g)$name[bn_nodes(g)$type == "gene"]
  n_ge <- sum(bn_edges(g)$from %in% genes)
  for (m in c("bernoulli", "pathway_size")) {
    r <- bn_rewire(g, m, seed = 3)
    rg <- bn_nodes(r)$name[bn_nodes(r)$type == "gene"]
    expect_equal(sum(bn_edges(r)$from %in% rg), n_ge)
    expect_equal(attr(r, "null_model"), m)
  }
  ps <- bn_rewire(g, "pathway_size", seed = 3)
  ge_size <- function(gr) {
    e <- bn_edges(gr); gn <- bn_nodes(gr)$name[bn_nodes(gr)$type == "gene"]
    table(e$to[e$from %in% gn])
  }
  pw <- names(ge_size(g))
  expect_equal(as.vector(ge_size(ps)[pw]), as.vector(ge_size(g)[pw]))
})

test_that("dataset helpers read from the cache when present", {
  skip_if_not(file.exists(file.path(binnr_cache_dir(), "hartman_aki", "sample_datamatrix.csv")),
              "Hartman AKI data not cached")
  aki <- hartman_aki()
  expect_equal(dim(aki$proteomics), c(197L, 554L))
  expect_false(anyNA(aki$proteomics))
  expect_equal(levels(aki$design$group), c("1", "2"))
  aki_na <- hartman_aki(impute = "none")
  expect_true(anyNA(aki_na$proteomics))
  aki_min <- hartman_aki(impute = "min")
  expect_true(all(aki_min$proteomics[is.na(aki_na$proteomics)] > 0))
})
