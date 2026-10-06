# binnr vs the Python `binn` package on Hartman et al.'s septic-AKI proteomics.
# Same data, same 5x5 folds (exported from Python), three binnr models:
#   A. "binn-exact": binn's own layered graph (4 strata, copy nodes), masked
#      dense, tanh, dropout 0.2 -- the stratified construction, fixed epochs
#   B. same graph, binnr defaults (early stopping on a validation split)
#   C. "topological": binnr's native DAG from Reactome (UniProt), full depth
# plus a lasso baseline. Then explanations: binn's SHAP vs binnr's
# gradient x activation importance on the same training folds.
suppressPackageStartupMessages({
  library(binnr)
  library(dplyr); library(tidyr); library(glmnet); library(yardstick)
})
set.seed(1)
out <- "Sys.getenv("BINNR_REPLICATION_OUT", "replication-out")"
aki <- hartman_aki(impute = "zero")           # as binn's dataloader does
X <- aki$proteomics; y <- aki$design$group    # levels "1","2"; "2" is binn's class 1
folds <- read.csv(file.path(out, "binn_pred.csv")) |> distinct(rep, fold, sample)
pred_binn <- read.csv(file.path(out, "binn_pred.csv"))

# ---- Graph A/B: binn's layered graph, copies made explicit --------------------
e <- read.csv(file.path(out, "binn_edges.csv")) |> filter(layer <= 3)
inputs <- read.csv(file.path(out, "binn_inputs.csv"))$protein
e <- e |> mutate(from_n = ifelse(layer == 0, from, paste0(from, "@L", layer)),
                 to_n = paste0(to, "@L", layer + 1))
g_binn <- bn_graph(data.frame(from = e$from_n, to = e$to_n), genes = inputs)
# a pathway appearing in more than one stratum: later appearances are copies
base <- sub("@L[0-9]+$", "", bn_nodes(g_binn)$name)
lvl <- as.integer(sub("^.*@L", "", bn_nodes(g_binn)$name))
first <- tapply(lvl, base, min)
dummies <- bn_nodes(g_binn)$name[bn_nodes(g_binn)$type == "pathway" & lvl > first[base]]
g_binn <- binnr:::mark_dummies(g_binn, dummies)
print(g_binn); cat("copy nodes:", length(dummies), "\n")

# ---- Graph C: binnr topological DAG ------------------------------------------------
g_full <- reactome_graph(identifiers = "uniprot")
g_topo <- bn_prune(g_full, genes = colnames(X))
print(g_topo)
genes_topo <- bn_nodes(g_topo)$name[bn_nodes(g_topo)$type == "gene"]
cat("proteins: binn graph", length(inputs), "| topological", length(genes_topo),
    "| shared", length(intersect(inputs, genes_topo)), "\n")

arch_binn <- bn_arch("binn", dropout = 0.2)
# binn trains 100 epochs, Adam 1e-4, batch 32 *with BatchNorm*; without BatchNorm binnr
# needs lr 1e-3 to reach a comparable training loss in 100 epochs (checked: 1e-4 underfits)
ctrl_fixed <- function(seed) bn_control(epochs = 100L, lr = 1e-3, batch_size = 32L, weight_decay = 0, validation = 0, patience = 100L, seed = seed)
ctrl_default <- function(seed) bn_control(seed = seed)

score <- function(method, rep, fold, prob, truth, secs) {
  tibble(method = method, rep = rep, fold = fold,
         auroc = roc_auc_vec(truth, prob, event_level = "second"),
         accuracy = accuracy_vec(truth, factor(ifelse(prob > 0.5, "2", "1"), levels(truth))),
         seconds = secs)
}
res <- list(); imp <- list()
for (r in 1:5) for (f in 1:5) {
  te <- which(aki$design$sample %in% folds$sample[folds$rep == r & folds$fold == f])
  tr <- setdiff(seq_along(y), te)
  seed <- r * 100 + f
  message(sprintf("rep %d fold %d", r, f))
  # Python binn, same fold
  pb <- pred_binn |> filter(rep == r, fold == f)
  pb <- pb[match(aki$design$sample[te], pb$sample), ]
  res[[length(res) + 1]] <- score("binn (Python, SHAP paper code)", r, f, pb$prob, y[te], NA)
  # A. binn-exact graph, fixed 100 epochs
  t0 <- Sys.time()
  fitA <- suppressMessages(bn_fit(X[tr, inputs], y[tr], graph = g_binn, arch = arch_binn, control = ctrl_fixed(seed)))
  pA <- predict(fitA, X[te, inputs], type = "prob")$.pred_2
  res[[length(res) + 1]] <- score("binnr, binn graph, fixed epochs", r, f, pA, y[te], as.numeric(Sys.time() - t0, units = "secs"))
  # B. binn-exact graph, binnr defaults
  t0 <- Sys.time()
  fitB <- suppressMessages(bn_fit(X[tr, inputs], y[tr], graph = g_binn, arch = arch_binn, control = ctrl_default(seed)))
  pB <- predict(fitB, X[te, inputs], type = "prob")$.pred_2
  res[[length(res) + 1]] <- score("binnr, binn graph, early stopping", r, f, pB, y[te], as.numeric(Sys.time() - t0, units = "secs"))
  # C. topological graph, binnr defaults
  t0 <- Sys.time()
  fitC <- suppressMessages(bn_fit(X[tr, ], y[tr], graph = g_topo, arch = arch_binn, control = ctrl_default(seed)))
  pC <- predict(fitC, X[te, ], type = "prob")$.pred_2
  res[[length(res) + 1]] <- score("binnr, topological graph", r, f, pC, y[te], as.numeric(Sys.time() - t0, units = "secs"))
  # lasso on the same proteins
  t0 <- Sys.time()
  m <- cv.glmnet(X[tr, ], y[tr], family = "binomial", nfolds = 5)
  pL <- as.numeric(predict(m, X[te, ], s = "lambda.min", type = "response"))
  res[[length(res) + 1]] <- score("Lasso (glmnet)", r, f, pL, y[te], as.numeric(Sys.time() - t0, units = "secs"))
  if (r == 1) {
    imp[[length(imp) + 1]] <- bind_rows(
      node_importance(fitA, X[tr, inputs]) |> mutate(model = "binnr, binn graph, fixed epochs"),
      node_importance(fitB, X[tr, inputs]) |> mutate(model = "binnr, binn graph, early stopping"),
      node_importance(fitC, X[tr, ]) |> mutate(model = "binnr, topological graph")
    ) |> mutate(rep = r, fold = f)
  }
}
res <- bind_rows(res); imp <- bind_rows(imp)
saveRDS(list(results = res, importance = imp, g_binn = g_binn, g_topo = g_topo), file.path(out, "binnr_compare.rds"))
print(res |> group_by(method) |> summarise(across(c(auroc, accuracy), list(mean = mean, sd = sd)), .groups = "drop"), n = 10, width = 120)
