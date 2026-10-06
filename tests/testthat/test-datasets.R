# Dataset helpers and downloaded knowledge graphs. Everything here reads from
# the cache and is skipped when the files are absent (no network in tests).
cache_has <- function(...) file.exists(file.path(binnr_cache_dir(), ...))

test_that("binnr_cache_dir honours BINNR_CACHE", {
  withr::with_envvar(c(BINNR_CACHE = "/tmp/some-cache"), expect_equal(binnr_cache_dir(), "/tmp/some-cache"))
  withr::with_envvar(c(BINNR_CACHE = ""), expect_equal(binnr_cache_dir(), tools::R_user_dir("binnr", "cache")))
  expect_identical(tcga_cache_dir(), binnr_cache_dir())
})

test_that("tcga_cohort validates input and reads the cache", {
  expect_error(tcga_cohort("LUAD"), "BRCA")
  skip_if_not(tcga_cached("PRAD"), "TCGA PRAD cohort not cached")
  prad <- tcga_cohort("PRAD")
  expect_named(prad, c("omics", "clinical"))
  expect_equal(dim(prad$omics$rna), c(491L, 1000L))
  expect_equal(names(prad$omics), c("rna", "cnv"))
  expect_setequal(colnames(prad$omics$rna), colnames(prad$omics$cnv))
  expect_equal(nrow(prad$clinical), 491L)
  expect_equal(levels(prad$clinical$grade), c("low", "high"))
  expect_match(attr(prad, "source"), "curatedTCGAData")
})

test_that("reactome_graph builds symbol and UniProt graphs from cached files", {
  skip_if_not(cache_has("reactome", "ReactomePathways.gmt.zip") && cache_has("reactome", "UniProt2Reactome.txt"),
              "Reactome files not cached")
  g <- reactome_graph(min_genes = 5)
  expect_s3_class(g, "bn_graph")
  n <- bn_nodes(g)
  expect_true(all(startsWith(n$name[n$type == "pathway"], "R-HSA-")))
  expect_true("TP53" %in% n$name)
  expect_false(any(duplicated(bn_edges(g))))
  gu <- reactome_graph(identifiers = "uniprot")
  expect_true("P04637" %in% bn_nodes(gu)$name) # TP53
  expect_false(any(grepl("^R-MMU", bn_nodes(gu)$name)))
  expect_s3_class(igraph::as.igraph(g), "igraph")
  expect_false(inherits(igraph::as.igraph(g), "bn_graph"))
  expect_error(reactome_graph(identifiers = "ensembl"))
})

test_that("pnet_prostate reproduces P-NET's inputs from the cached archive", {
  skip_if_not(cache_has("pnet", "_database", "prostate", "processed", "response_paper.csv"), "P-NET archive not cached")
  pc <- pnet_prostate()
  expect_equal(names(pc$omics), c("mut", "cna_del", "cna_amp"))
  expect_equal(dim(pc$omics$mut), c(1011L, 9229L))
  expect_true(all(unlist(lapply(pc$omics, function(m) all(m %in% 0:1)))))
  expect_equal(as.integer(table(pc$clinical$response)), c(678L, 333L))
  expect_setequal(unique(pc$splits$set), c("train", "validation", "test"))
  expect_equal(sum(pc$splits$set == "test"), 102L)
  pg <- pnet_prostate(cna = "gistic", genes = c("TP53", "AR", "PTEN"))
  expect_equal(names(pg$omics), c("mut", "cna"))
  expect_equal(colnames(pg$omics$mut), c("AR", "PTEN", "TP53"))
  expect_true(all(pg$omics$cna %in% -2:2))
})
