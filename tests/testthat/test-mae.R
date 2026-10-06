test_that("MultiAssayExperiment data are extracted and fitted", {
  brca <- brca_or_skip()
  skip_if_no_torch()
  skip_if_not_installed("MultiAssayExperiment")
  i <- 1:120
  cd <- as.data.frame(brca$clinical[i, ])
  rownames(cd) <- cd$patient
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rna = t(brca$omics$rna[i, 1:100]),
                       cnv = t(brca$omics$cnv[i[-(1:10)], 1:100])),
    colData = S4Vectors::DataFrame(cd)
  )
  d <- mae_to_omics(mae)
  expect_named(d$omics, c("rna", "cnv"))
  expect_equal(nrow(d$omics$rna), 110)   # patients in both assays
  expect_equal(d$clinical$patient, rownames(d$omics$cnv))
  g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy)
  fit <- suppressMessages(bn_fit(age ~ 1, data = mae, graph = g,
                                 control = bn_control(epochs = 2)))
  expect_equal(fit$family, "regression")
  expect_equal(fit$network$modalities, c("rna", "cnv"))
  expect_error(bn_fit(mae, 1, graph = g), "formula")
})

test_that("bn_arch and bn_control validate their arguments", {
  expect_error(bn_arch(hidden_dim = 0), "positive")
  expect_error(bn_arch(dropout = 1), "dropout")
  expect_error(bn_arch(weights = "foo"))
  expect_error(bn_control(validation = 1), "validation")
  a <- bn_arch("gnn", hidden_dim = 2)
  expect_equal(a$hidden_dim, 2L)
  expect_equal(a$weights, "shared")
  expect_output(print(a), "gnn preset")
  expect_output(print(bn_control()), "AdamW")
})
