# Random-graph controls on two of Caranzano et al.'s benchmarks, in binnr.
# 1. Hartman AKI proteomics: 3x5-fold CV, real Reactome vs four nulls.
# 2. P-NET prostate (primary vs metastatic): P-NET's own train/validation/test
#    split, 5 seeds, real vs nulls; plus a lasso baseline.
suppressPackageStartupMessages({
  library(binnr)
  library(dplyr); library(yardstick); library(glmnet)
})
out <- "Sys.getenv("BINNR_REPLICATION_OUT", "replication-out")"
nulls <- c("degree", "pathway_size", "bernoulli", "genes")
auc <- function(truth, prob) roc_auc_vec(truth, prob, event_level = "second")

# ---- 1. Hartman -----------------------------------------------------------------
aki <- hartman_aki(); X <- aki$proteomics; y <- aki$design$group
g <- bn_prune(reactome_graph(identifiers = "uniprot"), genes = colnames(X))
folds <- read.csv(file.path(out, "binn_pred.csv")) |> distinct(rep, fold, sample)
res1 <- list()
for (r in 1:3) for (f in 1:5) {
  te <- which(aki$design$sample %in% folds$sample[folds$rep == r & folds$fold == f]); tr <- setdiff(seq_along(y), te)
  seed <- r * 100 + f
  graphs <- c(list(Reactome = g), setNames(lapply(nulls, function(m) bn_rewire(g, m, seed = seed)), nulls))
  for (nm in names(graphs)) {
    t0 <- Sys.time()
    fit <- suppressMessages(bn_fit(X[tr, ], y[tr], graph = graphs[[nm]], arch = bn_arch("binn", dropout = 0.2),
                                   control = bn_control(seed = seed)))
    p <- predict(fit, X[te, ], type = "prob")$.pred_2
    res1[[length(res1) + 1]] <- tibble(task = "Hartman AKI", graph = nm, rep = r, fold = f, auroc = auc(y[te], p),
                                       seconds = as.numeric(Sys.time() - t0, units = "secs"))
  }
  m <- cv.glmnet(X[tr, ], y[tr], family = "binomial", nfolds = 5)
  res1[[length(res1) + 1]] <- tibble(task = "Hartman AKI", graph = "lasso", rep = r, fold = f,
                                     auroc = auc(y[te], as.numeric(predict(m, X[te, ], s = "lambda.min", type = "response"))), seconds = NA)
  message("Hartman rep ", r, " fold ", f)
}
res1 <- bind_rows(res1); saveRDS(res1, file.path(out, "nulls_hartman.rds"))
print(res1 |> group_by(graph) |> summarise(auroc = mean(auroc), sd = sd(auroc), secs = mean(seconds)))

# ---- 2. P-NET -------------------------------------------------------------------
pc <- pnet_prostate()
gp <- bn_prune(reactome_graph(), genes = colnames(pc$omics$mut))
print(gp)
sp <- pc$splits; y <- pc$clinical$response; names(y) <- pc$clinical$sample
tr <- sp$sample[sp$set %in% c("train", "validation")]; te <- sp$sample[sp$set == "test"]
rows <- function(i) lapply(pc$omics, function(m) m[i, , drop = FALSE])
res2 <- list()
for (seed in 1:5) {
  graphs <- c(list(Reactome = gp), setNames(lapply(nulls, function(m) bn_rewire(gp, m, seed = seed)), nulls))
  for (nm in names(graphs)) {
    t0 <- Sys.time()
    fit <- suppressMessages(bn_fit(rows(tr), y[tr], graph = graphs[[nm]], arch = bn_arch("binn"),
                                   control = bn_control(seed = seed)))
    p <- predict(fit, rows(te), type = "prob")$.pred_metastatic
    res2[[length(res2) + 1]] <- tibble(task = "P-NET prostate", graph = nm, seed = seed, auroc = auc(y[te], p),
                                       aupr = pr_auc_vec(y[te], p, event_level = "second"),
                                       seconds = as.numeric(Sys.time() - t0, units = "secs"))
    message("P-NET seed ", seed, " ", nm, ": AUROC ", round(res2[[length(res2)]]$auroc, 3), " (", round(res2[[length(res2)]]$seconds), "s)")
    saveRDS(bind_rows(res2), file.path(out, "nulls_pnet.rds"))
  }
}
Xl <- do.call(cbind, pc$omics)
m <- cv.glmnet(Xl[tr, ], y[tr], family = "binomial", nfolds = 5)
pl <- as.numeric(predict(m, Xl[te, ], s = "lambda.min", type = "response"))
res2[[length(res2) + 1]] <- tibble(task = "P-NET prostate", graph = "lasso", seed = NA, auroc = auc(y[te], pl),
                                   aupr = pr_auc_vec(y[te], pl, event_level = "second"), seconds = NA)
res2 <- bind_rows(res2); saveRDS(res2, file.path(out, "nulls_pnet.rds"))
print(res2 |> group_by(graph) |> summarise(auroc = mean(auroc), sd = sd(auroc), aupr = mean(aupr), secs = mean(seconds)))
