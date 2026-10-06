g_brca <- function() {
  pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
                labels = reactome_tcga$pathways)
}

test_that("regression recovers signal on a toy graph", {
  skip_if_no_torch()
  d <- toy_data(n = 200)
  fit <- quiet_fit(d$x, d$y, graph = toy_layered(),
                   arch = bn_arch(dropout = 0),
                   control = bn_control(epochs = 150, lr = 0.02, seed = 1))
  expect_s3_class(fit, "bn_model")
  expect_equal(fit$family, "regression")
  p <- predict(fit, d$x)
  expect_named(p, ".pred")
  expect_gt(stats::cor(p$.pred, d$y), 0.8)
})

test_that("classification with multi-omics and covariates (formula interface)", {
  brca <- brca_or_skip()
  skip_if_no_torch()
  i <- which(stats::complete.cases(brca$clinical))[1:200]
  om <- lapply(brca$omics, function(m) m[i, 1:200])
  cl <- brca$clinical[i, ]
  fit <- quiet_fit(pam50 ~ age + stage, data = cl, omics = om, graph = g_brca(),
                   control = bn_control(epochs = 5, seed = 1))
  expect_equal(fit$family, "classification")
  expect_true(fit$network$gene_layer)
  expect_equal(fit$network$modalities, c("rna", "cnv", "mut"))
  expect_gt(fit$n_cov, 1)
  pc <- predict(fit, om, covariates = cl)
  expect_s3_class(pc$.pred_class, "factor")
  expect_equal(levels(pc$.pred_class), levels(cl$pam50))
  pp <- predict(fit, om, covariates = cl, type = "prob")
  expect_equal(unname(rowSums(as.matrix(pp))), rep(1, nrow(pp)), tolerance = 1e-5)
  expect_error(predict(fit, om), "covariates")
  expect_error(predict(fit, om["rna"], covariates = cl), "lack")
})

test_that("Cox survival models train and predict", {
  brca <- brca_or_skip()
  skip_if_no_torch()
  skip_if_not_installed("survival")
  i <- which(stats::complete.cases(brca$clinical[c("age", "os_time")]))[1:300]
  om <- list(rna = brca$omics$rna[i, 1:300])
  cl <- brca$clinical[i, ]
  fit <- suppressMessages(bn_fit(survival::Surv(os_time, os_event) ~ age, data = cl,
                                 omics = om, graph = g_brca(),
                                 arch = bn_arch("gnn", hidden_dim = 2),
                                 control = bn_control(epochs = 5, seed = 1)))
  expect_equal(fit$family, "cox")
  expect_null(fit$params$head.bias)
  lp <- predict(fit, om, covariates = cl)
  expect_named(lp, ".pred_linear_pred")
  rk <- predict(fit, om, covariates = cl, type = "risk")
  expect_equal(rk$.pred_risk, exp(lp$.pred_linear_pred), tolerance = 1e-6)
  s <- predict(fit, om, covariates = cl, type = "survival", eval_time = c(365, 1825))
  surv1 <- s$.pred[[1]]
  expect_equal(surv1$.eval_time, c(365, 1825))
  expect_true(all(diff(surv1$.pred_survival) <= 0))
  expect_error(predict(fit, om, covariates = cl, type = "survival"), "eval_time")
})

test_that("cox partial likelihood matches survival::coxph", {
  skip_if_no_torch()
  skip_if_not_installed("survival")
  set.seed(1)
  n <- 50
  eta <- rnorm(n)
  time <- rexp(n) + 0.01
  status <- rbinom(n, 1, 0.7)
  # coxph log-likelihood at fixed coefficients via offset
  cf <- survival::coxph(survival::Surv(time, status) ~ offset(eta), ties = "breslow")
  ours <- binnr:::cox_partial_loss(torch::torch_tensor(eta), time,
                                      torch::torch_tensor(as.numeric(status)))$item()
  expect_equal(-ours * sum(status), cf$loglik[1], tolerance = 1e-4)
})

test_that("fits are reproducible with a seed and survive saveRDS", {
  skip_if_no_torch()
  d <- toy_data(n = 60)
  ctrl <- bn_control(epochs = 5, seed = 11)
  f1 <- quiet_fit(d$x, d$y, graph = toy_ragged(), control = ctrl)
  f2 <- quiet_fit(d$x, d$y, graph = toy_ragged(), control = ctrl)
  expect_equal(predict(f1, d$x), predict(f2, d$x))
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(f1, path)
  f3 <- readRDS(path)
  expect_equal(predict(f3, d$x), predict(f1, d$x))
})

test_that("predict does not disturb the user's RNG", {
  skip_if_no_torch()
  d <- toy_data(n = 30)
  f <- quiet_fit(d$x, d$y, graph = toy_layered(), control = bn_control(epochs = 2))
  set.seed(5); a <- runif(1)
  set.seed(5); invisible(predict(f, d$x)); b <- runif(1)
  expect_equal(a, b)
})

test_that("informative errors for bad inputs", {
  skip_if_no_torch()
  d <- toy_data(n = 20)
  g <- toy_layered()
  expect_error(bn_fit(d$x, d$y[-1], graph = g), "has 19")
  x_na <- d$x; x_na[1, 1] <- NA
  expect_error(bn_fit(x_na, d$y, graph = g), "missing")
  expect_error(bn_fit(d$x, d$y, graph = list()), "bn_graph")
  x_bad <- d$x; colnames(x_bad) <- paste0("X", 1:4)
  expect_error(suppressMessages(bn_fit(x_bad, d$y, graph = g)), "None")
  expect_error(bn_fit(d$x, as.Date(1:20), graph = g), "outcome")
  expect_error(bn_fit(list(d$x), d$y, graph = g), "names")
  expect_error(predict(quiet_fit(d$x, d$y, graph = g, control = bn_control(epochs = 1))), "new_data")
})

test_that("unmapped features are dropped with a message", {
  skip_if_no_torch()
  d <- toy_data(n = 30, genes = c("G1", "G2", "G3", "G4", "XX"))
  expect_message(f <- bn_fit(d$x, d$y, graph = toy_layered(),
                             control = bn_control(epochs = 1)),
                 "4 of 5")
  expect_equal(f$network$n_input, 4)
})

test_that("methods: print, glance, tidy, autoplot", {
  skip_if_no_torch()
  d <- toy_data(n = 40)
  f <- quiet_fit(d$x, d$y, graph = toy_ragged(), control = bn_control(epochs = 5, seed = 1))
  expect_output(print(f), "bn_model")
  gl <- glance(f)
  expect_equal(nrow(gl), 1)
  expect_equal(gl$n_features, 4)
  td <- tidy(f)
  expect_equal(nrow(td), length(f$network$src))
  expect_false(anyNA(td$weight))
  skip_if_not_installed("ggplot2")
  expect_s3_class(ggplot2::autoplot(f), "ggplot")
})

test_that("node importance ranks informative genes highest", {
  skip_if_no_torch()
  d <- toy_data(n = 200)
  f <- quiet_fit(d$x, d$y, graph = toy_ragged(), arch = bn_arch(dropout = 0),
                 control = bn_control(epochs = 100, lr = 0.02, seed = 1))
  imp <- node_importance(f, d$x)
  expect_s3_class(imp, "bn_importance")
  expect_equal(nrow(imp), nrow(f$network$nodes))
  inp <- imp[imp$type == "input", ]
  expect_setequal(inp$node[1:2], c("G1|x", "G3|x"))
  skip_if_not_installed("ggplot2")
  expect_s3_class(ggplot2::autoplot(imp, node_type = "input"), "ggplot")
})

test_that("gene layer can be switched off and multi-dimensional embeddings work", {
  skip_if_no_torch()
  d <- toy_data(n = 40)
  om <- list(a = d$x, b = d$x + 1)
  f <- quiet_fit(om, d$y, graph = toy_layered(),
                 arch = bn_arch(gene_layer = FALSE, hidden_dim = 3, readout = "pathways"),
                 control = bn_control(epochs = 2))
  expect_false(any(f$network$nodes$type == "gene"))
  expect_equal(dim(f$params$embed), c(2L, 3L))
  expect_equal(nrow(predict(f, om)), 40)
})

test_that("Cox loss handles tied times like coxph(ties = 'breslow')", {
  skip_if_no_torch()
  skip_if_not_installed("survival")
  set.seed(3)
  n <- 60
  eta <- rnorm(n)
  time <- sample(1:8, n, replace = TRUE)
  status <- rbinom(n, 1, 0.6)
  cf <- survival::coxph(survival::Surv(time, status) ~ offset(eta), ties = "breslow")
  ours <- binnr:::cox_partial_loss(torch::torch_tensor(eta), time,
                                      torch::torch_tensor(as.numeric(status)))$item()
  expect_equal(-ours * sum(status), cf$loglik[1], tolerance = 1e-4)
})

test_that("single-modality models accept a matrix or a one-element list", {
  skip_if_no_torch()
  d <- toy_data(n = 30)
  f <- quiet_fit(list(rna = d$x), d$y, graph = toy_layered(), control = bn_control(epochs = 2))
  expect_equal(predict(f, d$x), predict(f, list(rna = d$x)))
  expect_equal(predict(f, list(other = d$x)), predict(f, d$x))
})

test_that("formula method ignores NAs in unmapped columns", {
  skip_if_no_torch()
  d <- toy_data(n = 40, genes = c("G1", "G2", "G3", "G4", "JUNK"))
  d$x[, "JUNK"] <- NA
  dat <- data.frame(y = d$y)
  f <- quiet_fit(y ~ 1, data = dat, omics = d$x, graph = toy_layered(),
                 control = bn_control(epochs = 2))
  expect_equal(f$n, 40)
  dat$y[1:35] <- NA
  expect_error(quiet_fit(y ~ 1, data = dat, omics = d$x, graph = toy_layered()), "at least 10")
})

test_that("tidy and glance are registered S3 methods", {
  expect_true(!is.null(utils::getS3method("tidy", "bn_model", optional = TRUE,
                                          envir = asNamespace("generics"))))
  expect_true(!is.null(utils::getS3method("glance", "bn_model", optional = TRUE,
                                          envir = asNamespace("generics"))))
})

test_that("covariate edge cases give clear messages", {
  skip_if_no_torch()
  d <- toy_data(n = 40)
  cov <- data.frame(`ER status` = rep(c("pos", "neg"), 20), check.names = FALSE)
  f <- quiet_fit(d$x, d$y, graph = toy_layered(), covariates = cov,
                 control = bn_control(epochs = 2))
  expect_equal(f$n_cov, 1)
  bad <- data.frame(`ER status` = rep("unknown", 40), check.names = FALSE)
  expect_error(predict(f, d$x, covariates = bad), "incompatible")
  f2 <- quiet_fit(d$x, d$y, graph = toy_layered(), control = bn_control(epochs = 2))
  expect_warning(predict(f2, d$x, covariates = cov), "no covariates")
})

test_that("predict leaves torch's RNG untouched", {
  skip_if_no_torch()
  d <- toy_data(n = 30)
  f <- quiet_fit(d$x, d$y, graph = toy_layered(), control = bn_control(epochs = 2))
  torch::torch_manual_seed(9); a <- torch::torch_rand(1)$item()
  torch::torch_manual_seed(9); invisible(predict(f, d$x)); b <- torch::torch_rand(1)$item()
  expect_equal(a, b)
})

test_that("survival outcomes without events are rejected", {
  skip_if_no_torch()
  skip_if_not_installed("survival")
  d <- toy_data(n = 20)
  expect_error(bn_fit(d$x, survival::Surv(1:20, rep(0, 20)), graph = toy_layered()), "no events")
})

test_that("update() refits with a different graph", {
  skip_if_no_torch()
  d <- toy_data(n = 30)
  x <- d$x
  y <- d$y
  f <- suppressMessages(bn_fit(x, y, graph = toy_layered(), control = bn_control(epochs = 2)))
  expect_equal(f$call[[1]], quote(bn_fit))
  f2 <- suppressMessages(update(f, graph = toy_ragged()))
  expect_s3_class(f2, "bn_model")
  expect_true("P4" %in% f2$network$nodes$name)
})

test_that("node_activations returns pathway scores in wide and long form", {
  skip_if_no_torch()
  d <- toy_data(n = 30)
  rownames(d$x) <- paste0("s", 1:30)
  f <- quiet_fit(d$x, d$y, graph = toy_ragged(), control = bn_control(epochs = 2))
  a <- node_activations(f, d$x, type = "pathway")
  expect_equal(dim(a), c(30L, 4L))
  expect_equal(rownames(a), paste0("s", 1:30))
  expect_true(all(abs(a) <= 1))  # tanh
  l <- node_activations(f, d$x, format = "long")
  expect_equal(nrow(l), 30 * 4)
  f3 <- quiet_fit(d$x, d$y, graph = toy_ragged(), arch = bn_arch(hidden_dim = 3),
                  control = bn_control(epochs = 2))
  expect_equal(ncol(node_activations(f3, d$x, type = "pathway")), 12L)
})

test_that("bn_tune selects the configuration with the lowest validation loss", {
  skip_if_no_torch()
  d <- toy_data(n = 60)
  fit <- bn_tune(d$x, d$y, graph = toy_ragged(),
                 grid = expand.grid(hidden_dim = c(1, 2), max_level = c(1, Inf)),
                 control = bn_control(epochs = 3, seed = 1))
  expect_s3_class(fit, "bn_model")
  expect_equal(nrow(fit$tuning), 4)
  expect_equal(sum(fit$tuning$selected), 1)
  expect_equal(fit$best_loss, min(fit$tuning$validation_loss))
  # max_level = 1 removes the level-2 pathways P3, P4
  expect_error(bn_tune(d$x, d$y, graph = toy_ragged(), grid = data.frame(foo = 1)), "Unknown")
})

test_that("baseline hazard matches survival::basehaz and survival predictions are valid", {
  skip_if_not_installed("survival")
  set.seed(3)
  lp <- rnorm(40); time <- rexp(40, exp(lp)); status <- rbinom(40, 1, 0.7)
  time[1:4] <- time[5] # ties
  bl <- binnr:::breslow(lp, time, status)
  cf <- survival::coxph(survival::Surv(time, status) ~ offset(lp), ties = "breslow")
  bh <- survival::basehaz(cf, centered = FALSE)
  expect_equal(bl$cumhaz, bh$hazard)
  expect_true(all(diff(bl$cumhaz) >= 0))
})

test_that("importance stability and targeted importance work on a toy graph", {
  skip_if_no_torch()
  d <- toy_data(n = 120)
  g <- toy_ragged()
  f <- quiet_fit(d$x, d$y, graph = g, control = bn_control(epochs = 15, seed = 1))
  st <- bn_importance_stability(f, d$x, d$y, graph = g, seeds = 1:2)
  expect_true(all(c("node", "mean_importance", "sd_importance", "mean_rank", "sd_rank") %in% names(st)))
  expect_setequal(unique(attr(st, "scores")$seed), 1:2)
  expect_equal(st$mean_rank, sort(st$mean_rank))
  imp <- node_importance(f, d$x, normalize = TRUE)
  expect_equal(max(imp$importance), 1)
  # classification: importance for one class
  yc <- factor(ifelse(d$y > 0, "hi", "lo"))
  fc <- quiet_fit(d$x, yc, graph = g, control = bn_control(epochs = 10, seed = 1))
  imp_hi <- node_importance(fc, d$x, target = "hi")
  expect_s3_class(imp_hi, "bn_importance")
  expect_error(node_importance(fc, d$x, target = "nope"), "target")
})

test_that("as_omics rejects malformed inputs with clear messages", {
  g <- toy_layered()
  expect_error(bn_fit(list(), 1:3, graph = g), "list")
  expect_error(bn_fit(list(a = matrix(1, 3, 2), a = matrix(1, 3, 2)), rnorm(3), graph = g), "unique")
  expect_error(bn_fit(list(rna = matrix("a", 3, 2)), rnorm(3), graph = g), "numeric")
  expect_error(bn_fit(list(rna = matrix(1, 3, 2)), rnorm(3), graph = g), "column names")
  expect_error(bn_fit(list(rna = matrix(1, 3, 2, dimnames = list(NULL, c("G1", "G2"))),
                           cnv = matrix(1, 4, 2, dimnames = list(NULL, c("G1", "G2")))), rnorm(3), graph = g),
               "same number of rows")
})
