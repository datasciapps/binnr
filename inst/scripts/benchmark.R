# Benchmark binnr against standard R learners on the example data.
#
# Three tasks, each evaluated by repeated stratified k-fold cross-validation
# with identical folds for every method:
#   * BRCA: PAM50 subtype (4 classes) from RNA + CNV + mutation
#   * BRCA: overall survival (Cox) from RNA + CNV + mutation + age + stage
#   * PRAD: Gleason grade >= 8 vs <= 7 (binary) from RNA + CNV, as in P-NET
#
# Hyperparameters are chosen *inside* each training fold with a matched
# budget: binnr models by validation loss over a small grid (graph depth,
# state dimension, weight decay) via bn_tune(); glmnet's penalty and xgboost's
# number of rounds by inner 5-fold CV; ranger at its defaults. Null,
# incomplete and padded graphs reuse the configuration selected for the
# Reactome BINN and are re-drawn for every split, so their randomness is part
# of the resampling variability.
#
# Usage (from the package root):
#   Rscript inst/scripts/benchmark.R [task] [folds] [repeats] [rep]
#     task    one of pam50, survival, gleason, all (default all)
#     folds   k (default 5);  repeats  number of repeats (default 5)
#     rep     optional: run only this repeat (for parallel jobs)
#   Rscript inst/scripts/benchmark.R merge
# Each (task, repeat) is written to inst/extdata/parts/<task>_rep<r>.rds and
# skipped if it already exists; `merge` combines the parts into
# inst/extdata/benchmark.rds, which the vignette and the abstract read.
# Set BINNR_THREADS to the number of CPU cores available to this job.

`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages({
  library(binnr)
  library(survival)
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(yardstick)
})
# TCGA example cohorts (downloaded once and cached; see ?tcga_cohort)
brca <- tcga_cohort("BRCA")
prad <- tcga_cohort("PRAD")

args <- commandArgs(trailingOnly = TRUE)
task <- if (length(args) >= 1) args[[1]] else "all"
k <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
reps <- if (length(args) >= 3) as.integer(args[[3]]) else 5L
only_rep <- if (length(args) >= 4) as.integer(args[[4]]) else NA_integer_
out_dir <- file.path("inst", "extdata")
parts_dir <- file.path(out_dir, "parts")
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)
n_threads <- as.integer(Sys.getenv("BINNR_THREADS", "2"))
torch::torch_set_num_threads(n_threads)

graph <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
                       labels = reactome_tcga$pathways)

# ---- Tuning grids ---------------------------------------------------------------
# Selected by validation loss inside the training fold (bn_tune()).
binn_grid <- expand.grid(max_level = c(5, Inf), weight_decay = c(1e-3, 1e-2))
gnn_grid <- expand.grid(max_level = c(5, Inf), hidden_dim = c(4L, 8L))

# Graph variants of the tuned BINN: each takes the cohort graph (pruned to the
# selected depth) and the split (r, f)
graph_variants <- list(
  "BINN, stratified (padded)" = function(g, r, f) bn_pad(g),
  "BINN, rewired graph" = function(g, r, f) bn_rewire(g, "degree", seed = 1000 * r + f),
  "BINN, random graph" = function(g, r, f) bn_rewire(g, "bernoulli", seed = 1000 * r + f),
  "BINN, 50% of annotations removed" = function(g, r, f) bn_drop_edges(g, 0.5, seed = 1000 * r + f)
)

# ---- Helpers ------------------------------------------------------------------
folds_for <- function(strata, k, rep) {
  withr::with_seed(1000 + rep, {
    f <- integer(length(strata))
    for (s in unique(strata)) {
      i <- which(strata == s)
      f[i] <- sample(rep_len(seq_len(k), length(i)))
    }
    f
  })
}

rows_of <- function(omics, i) lapply(omics, function(m) m[i, , drop = FALSE])
flatten <- function(omics, i) do.call(cbind, Map(function(m, nm) {
  x <- m[i, , drop = FALSE]
  colnames(x) <- paste(colnames(x), nm, sep = "_")
  x
}, omics, names(omics)))

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(value = val, seconds = proc.time()[["elapsed"]] - t0)
}

xgb_tuned <- function(params, dm, ...) {
  cvm <- xgb.cv(params = params, data = dm, nrounds = 600, nfold = 5,
                early_stopping_rounds = 30, verbose = 0)
  # xgboost >= 3 stores the result in $early_stop; older versions at top level
  best <- cvm$early_stop$best_iteration %||% cvm$best_iteration %||% cvm$niter
  xgb.train(params = params, data = dm, nrounds = max(best, 10L), verbose = 0)
}

# Fit the tuned BINN, its graph variants and the tuned GNN. `fit_fun(graph,
# arch, ctrl)`, `tune_fun(grid, arch, ctrl)` and `pred_fun(fit)` are supplied
# by the task; `score(method, pred, seconds)` returns a data.frame.
run_binnr <- function(g, r, f, fit_fun, tune_fun, pred_fun, score) {
  res <- list()
  seed <- r * 100 + f
  t <- timed({
    binn <- tune_fun(binn_grid, bn_arch("binn"), bn_control(seed = seed))
    pred_fun(binn)
  })
  res[[length(res) + 1]] <- score("BINN (binnr)", t$value, t$seconds)
  sel <- binn$tuning[binn$tuning$selected, ]
  g_sel <- if (is.finite(sel$max_level)) bn_prune(g, max_level = sel$max_level) else g
  ctrl_sel <- bn_control(seed = seed, weight_decay = sel$weight_decay)
  for (v in names(graph_variants)) {
    t <- timed(pred_fun(fit_fun(graph_variants[[v]](g_sel, r, f), bn_arch("binn"), ctrl_sel)))
    res[[length(res) + 1]] <- score(v, t$value, t$seconds)
  }
  t <- timed({
    gnn <- tune_fun(gnn_grid, bn_arch("gnn"), bn_control(seed = seed))
    pred_fun(gnn)
  })
  res[[length(res) + 1]] <- score("GNN (binnr)", t$value, t$seconds)
  # record the selected configurations
  selected <- rbind(
    data.frame(method = "BINN (binnr)", rep = r, fold = f, max_level = sel$max_level,
               hidden_dim = 1L, weight_decay = sel$weight_decay),
    with(gnn$tuning[gnn$tuning$selected, ],
         data.frame(method = "GNN (binnr)", rep = r, fold = f, max_level = max_level,
                    hidden_dim = hidden_dim, weight_decay = 1e-3))
  )
  for (i in seq_along(res)) attr(res[[i]], "selected") <- NULL
  attr(res, "selected") <- selected
  res
}

bind_results <- function(res) {
  out <- do.call(rbind, res)
  attr(out, "selected") <- do.call(rbind, lapply(res, attr, "selected"))
  out
}

# ---- Classification (PAM50, Gleason) ---------------------------------------------
run_classification <- function(task_name, omics, y, binary, rep_seq) {
  idx <- which(!is.na(y))
  y <- droplevels(y[idx])
  X <- flatten(omics, idx)
  g <- bn_prune(graph, genes = colnames(omics[[1]]))
  res <- list()
  selected <- list()
  for (r in rep_seq) {
    fold <- folds_for(y, k, r)
    for (f in seq_len(k)) {
      tr <- which(fold != f); te <- which(fold == f)
      message(sprintf("%s rep %d fold %d", task_name, r, f))
      # `prob` is an n_test x n_class matrix of class probabilities
      score <- function(method, prob, seconds) {
        prob <- prob[, levels(y), drop = FALSE]
        pred <- factor(levels(y)[max.col(prob, ties.method = "first")], levels = levels(y))
        truth <- y[te]
        # yardstick conventions: brier_class() is half the multiclass Brier
        # score; the second factor level is the event for roc_auc()
        metrics <- c(accuracy = yardstick::accuracy_vec(truth, pred),
                     balanced_accuracy = yardstick::bal_accuracy_vec(truth, pred),
                     brier = if (binary) {
                       yardstick::brier_class_vec(truth, prob[, 2], event_level = "second")
                     } else {
                       yardstick::brier_class_vec(truth, prob)
                     })
        if (binary) {
          metrics <- c(metrics, auroc = yardstick::roc_auc_vec(truth, prob[, 2], event_level = "second"))
        }
        data.frame(task = task_name, method = method, rep = r, fold = f,
                   metric = names(metrics), value = unname(metrics), seconds = seconds)
      }
      fam <- if (binary) "binomial" else "multinomial"
      t <- timed({
        m <- cv.glmnet(X[tr, ], y[tr], family = fam, alpha = 1, nfolds = 5)
        p <- predict(m, X[te, ], s = "lambda.min", type = "response")
        if (binary) cbind(1 - p[, 1], p[, 1], deparse.level = 0) else p[, , 1]
      })
      colnames(t$value) <- levels(y)
      res[[length(res) + 1]] <- score("Lasso (glmnet)", t$value, t$seconds)
      t <- timed({
        m <- ranger(x = X[tr, ], y = y[tr], num.trees = 500, num.threads = n_threads,
                    seed = r, probability = TRUE)
        predict(m, X[te, ])$predictions
      })
      res[[length(res) + 1]] <- score("Random forest (ranger)", t$value, t$seconds)
      t <- timed({
        d <- xgb.DMatrix(X[tr, ], label = as.integer(y[tr]) - 1L)
        params <- list(eta = 0.05, max_depth = 3, subsample = 0.8,
                       colsample_bytree = 0.3, nthread = n_threads)
        if (binary) {
          params$objective <- "binary:logistic"; params$eval_metric <- "logloss"
        } else {
          params$objective <- "multi:softprob"; params$num_class <- nlevels(y)
          params$eval_metric <- "mlogloss"
        }
        m <- xgb_tuned(params, d)
        p <- predict(m, X[te, ])
        if (binary) cbind(1 - p, p) else {
          if (is.null(dim(p))) p <- matrix(p, ncol = nlevels(y), byrow = TRUE)
          p
        }
      })
      colnames(t$value) <- levels(y)
      res[[length(res) + 1]] <- score("Boosting (xgboost)", t$value, t$seconds)
      fit_fun <- function(graph, arch, ctrl) {
        suppressMessages(bn_fit(rows_of(omics, idx[tr]), y[tr], graph = graph,
                                arch = arch, control = ctrl))
      }
      tune_fun <- function(grid, arch, ctrl) {
        bn_tune(rows_of(omics, idx[tr]), y[tr], graph = g, grid = grid,
                arch = arch, control = ctrl)
      }
      pred_fun <- function(fit) {
        p <- as.matrix(predict(fit, rows_of(omics, idx[te]), type = "prob"))
        colnames(p) <- sub("^\\.pred_", "", colnames(p))
        p
      }
      b <- run_binnr(g, r, f, fit_fun, tune_fun, pred_fun, score)
      selected[[length(selected) + 1]] <- cbind(task = task_name, attr(b, "selected"))
      res <- c(res, b)
    }
  }
  out <- do.call(rbind, res)
  attr(out, "selected") <- do.call(rbind, selected)
  out
}

# ---- Survival (BRCA) -----------------------------------------------------------
run_survival <- function(rep_seq) {
  cl <- brca$clinical
  omics <- brca$omics
  idx <- which(stats::complete.cases(cl[c("age", "stage", "os_time", "os_event")]))
  d <- as.data.frame(cl[idx, ])
  y <- Surv(d$os_time, d$os_event)
  Z <- stats::model.matrix(~ age + stage, d)[, -1]
  X <- cbind(flatten(omics, idx), Z)
  pf <- c(rep(1, ncol(X) - ncol(Z)), rep(0, ncol(Z)))
  g <- bn_prune(graph, genes = colnames(omics[[1]]))
  res <- list()
  selected <- list()
  for (r in rep_seq) {
    fold <- folds_for(d$os_event, k, r)
    for (f in seq_len(k)) {
      tr <- which(fold != f); te <- which(fold == f)
      message(sprintf("Survival rep %d fold %d", r, f))
      score <- function(method, risk, seconds) {
        ci <- concordance(y[te] ~ risk, reverse = TRUE)$concordance
        data.frame(task = "Overall survival", method = method, rep = r, fold = f,
                   metric = "c_index", value = ci, seconds = seconds)
      }
      t <- timed({
        m <- coxph(Surv(os_time, os_event) ~ age + stage, data = d[tr, ])
        predict(m, d[te, ])
      })
      res[[length(res) + 1]] <- score("Cox, clinical only", t$value, t$seconds)
      t <- timed({
        m <- cv.glmnet(X[tr, ], y[tr], family = "cox", alpha = 1, nfolds = 5,
                       penalty.factor = pf)
        as.numeric(predict(m, X[te, ], s = "lambda.min"))
      })
      res[[length(res) + 1]] <- score("Lasso (glmnet)", t$value, t$seconds)
      t <- timed({
        m <- ranger(x = X[tr, ], y = y[tr], num.trees = 500, num.threads = n_threads, seed = r)
        rowSums(predict(m, X[te, ])$chf)
      })
      res[[length(res) + 1]] <- score("Random forest (ranger)", t$value, t$seconds)
      t <- timed({
        lab <- ifelse(d$os_event[tr] == 1, d$os_time[tr], -d$os_time[tr])
        params <- list(objective = "survival:cox", eta = 0.03, max_depth = 2,
                       subsample = 0.8, colsample_bytree = 0.3, nthread = n_threads)
        m <- xgb_tuned(params, xgb.DMatrix(X[tr, ], label = lab))
        predict(m, X[te, ])
      })
      res[[length(res) + 1]] <- score("Boosting (xgboost)", t$value, t$seconds)
      cov_tr <- d[tr, c("age", "stage")]
      cov_te <- d[te, c("age", "stage")]
      fit_fun <- function(graph, arch, ctrl) {
        suppressMessages(bn_fit(rows_of(omics, idx[tr]), y[tr], graph = graph,
                                covariates = cov_tr, arch = arch, control = ctrl))
      }
      tune_fun <- function(grid, arch, ctrl) {
        bn_tune(rows_of(omics, idx[tr]), y[tr], graph = g, grid = grid,
                covariates = cov_tr, arch = arch, control = ctrl)
      }
      pred_fun <- function(fit) {
        predict(fit, rows_of(omics, idx[te]), covariates = cov_te)$.pred_linear_pred
      }
      b <- run_binnr(g, r, f, fit_fun, tune_fun, pred_fun, score)
      selected[[length(selected) + 1]] <- cbind(task = "Overall survival", attr(b, "selected"))
      res <- c(res, b)
    }
  }
  out <- do.call(rbind, res)
  attr(out, "selected") <- do.call(rbind, selected)
  out
}

# ---- Run ----------------------------------------------------------------------
tasks <- c(pam50 = "PAM50 subtype", survival = "Overall survival", gleason = "Gleason grade")
run_task <- function(t, rep_seq) switch(t,
  pam50 = run_classification(tasks[["pam50"]], brca$omics, brca$clinical$pam50,
                             binary = FALSE, rep_seq = rep_seq),
  survival = run_survival(rep_seq),
  gleason = run_classification(tasks[["gleason"]], prad$omics, prad$clinical$grade,
                               binary = TRUE, rep_seq = rep_seq)
)
part_file <- function(t, r) file.path(parts_dir, sprintf("%s_rep%d.rds", t, r))

if (task == "merge") {
  files <- list.files(parts_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No parts found in ", parts_dir)
  parts <- lapply(files, readRDS)
  bench <- do.call(rbind, parts)
  attr(bench, "selected") <- do.call(rbind, lapply(parts, attr, "selected"))
  attr(bench, "info") <- list(
    folds = k, repeats = length(unique(bench$rep)), date = Sys.Date(),
    binnr = as.character(utils::packageVersion("binnr")),
    session = utils::sessionInfo()$R.version$version.string,
    tuning = list(binn = binn_grid, gnn = gnn_grid)
  )
  saveRDS(bench, file.path(out_dir, "benchmark.rds"))
  message("Merged ", length(files), " parts (", attr(bench, "info")$repeats,
          " repeats) into ", file.path(out_dir, "benchmark.rds"))
} else {
  todo <- if (task == "all") names(tasks) else task
  rep_seq <- if (is.na(only_rep)) seq_len(reps) else only_rep
  for (t in todo) for (r in rep_seq) {
    if (file.exists(part_file(t, r))) { message("Skipping ", t, " rep ", r); next }
    saveRDS(run_task(t, r), part_file(t, r))
  }
  message("Done. Merge with: Rscript inst/scripts/benchmark.R merge")
}
