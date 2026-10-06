# Compare explanations: binn's SHAP node importance vs binnr's |gradient x activation|
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
out <- "Sys.getenv("BINNR_REPLICATION_OUT", "replication-out")"
cmp <- readRDS(file.path(out, "binnr_compare.rds"))
res <- cmp$results; imp <- cmp$importance
# ---- performance table with QV comparison intervals and paired differences
k <- 5; J <- 25
qv <- function(d) {
  d$split <- interaction(d$rep, d$fold); d$method <- factor(d$method)
  fit <- lm(auroc ~ method + split, d); idx <- grep("^method", names(coef(fit)))
  m <- nlevels(d$method); V <- matrix(0, m, m); V[-1, -1] <- vcov(fit)[idx, idx]; V <- V * (1 + J / (k - 1))
  q <- qvcalc::qvcalc(V, estimates = c(0, coef(fit)[idx]), labels = levels(d$method))$qvframe
  tibble(method = levels(d$method), mean = as.numeric(tapply(d$auroc, d$method, mean)), qse = sqrt(q$quasiVar)) |>
    mutate(lo = mean - 1.39 * qse, hi = mean + 1.39 * qse)
}
perf <- qv(res) |> arrange(desc(mean))
acc <- res |> group_by(method) |> summarise(accuracy = mean(accuracy), seconds = mean(seconds, na.rm = TRUE))
print(left_join(perf, acc, by = "method"), width = 120)
paired <- function(a, b) {
  x <- res |> filter(method == a) |> select(rep, fold, a = auroc)
  y <- res |> filter(method == b) |> select(rep, fold, b = auroc)
  d <- inner_join(x, y, by = c("rep", "fold")) |> mutate(d = a - b)
  se <- sqrt((1 / J + 1 / (k - 1)) * var(d$d))
  sprintf("%s minus %s: %+.3f [%+.3f, %+.3f]", a, b, mean(d$d), mean(d$d) - qt(.975, J - 1) * se, mean(d$d) + qt(.975, J - 1) * se)
}
cat(paired("binnr, binn graph, fixed epochs", "binn (Python, SHAP paper code)"), "\n")
cat(paired("binnr, binn graph, early stopping", "binn (Python, SHAP paper code)"), "\n")
cat(paired("binnr, topological graph", "binnr, binn graph, early stopping"), "\n")
cat(paired("Lasso (glmnet)", "binn (Python, SHAP paper code)"), "\n")

# ---- explanations: node-level mean |SHAP| (binn) per fold of rep 1
shap <- read.csv(file.path(out, "binn_shap.csv"))
shap_node <- shap |> group_by(fold, source_layer, source_node) |>
  summarise(shap = mean(importance), .groups = "drop") |>
  mutate(base = sub("@L[0-9]+$", "", source_node))
# binnr: node importance per model/fold; collapse copy nodes to base pathway (max over copies)
imp_node <- imp |> mutate(type = ifelse(type == "input", "gene", type), node = sub("\\|x$", "", node), base = sub("@L[0-9]+$", "", node)) |>
  group_by(model, fold, base, type) |> summarise(importance = max(importance), .groups = "drop")

compare_rank <- function(model_name, layer_type) {
  d <- inner_join(
    imp_node |> filter(model == model_name, type %in% layer_type),
    shap_node |> group_by(fold, base) |> summarise(shap = max(shap), .groups = "drop"),
    by = c("fold", "base"))
  d |> group_by(fold) |> summarise(
    n = n(),
    spearman = cor(importance, shap, method = "spearman"),
    top10_overlap = length(intersect(base[order(-importance)][1:10], base[order(-shap)][1:10])),
    top20_overlap = length(intersect(base[order(-importance)][1:20], base[order(-shap)][1:20])),
    .groups = "drop") |> mutate(model = model_name, nodes = paste(layer_type, collapse = "/"))
}
agree <- bind_rows(
  compare_rank("binnr, binn graph, fixed epochs", "gene"),
  compare_rank("binnr, binn graph, early stopping", "gene"),
  compare_rank("binnr, topological graph", "gene"),
  compare_rank("binnr, binn graph, fixed epochs", c("pathway", "dummy")),
  compare_rank("binnr, binn graph, early stopping", c("pathway", "dummy")),
  compare_rank("binnr, topological graph", c("pathway", "dummy"))
)
cat("\nAgreement binnr importance vs binn SHAP (per fold, rep 1):\n")
print(agree |> group_by(model, nodes) |> summarise(n = mean(n), spearman = mean(spearman), top10 = mean(top10_overlap), top20 = mean(top20_overlap), .groups = "drop"), width = 120)

# stability across folds (within method): Spearman between folds of rank vectors
stab <- function(d, val) {
  w <- d |> group_by(fold, base) |> summarise(v = max(.data[[val]]), .groups = "drop") |>
    pivot_wider(names_from = fold, values_from = v) |> select(-base)
  cm <- cor(as.matrix(w), use = "pairwise", method = "spearman"); mean(cm[upper.tri(cm)])
}
cat("\nStability across the 5 folds (mean pairwise Spearman of node rankings):\n")
cat("binn SHAP, proteins:", round(stab(shap_node |> filter(source_layer == 0), "shap"), 2),
    " pathways:", round(stab(shap_node |> filter(source_layer > 0) |> group_by(fold, base) |> summarise(shap = max(shap), .groups = "drop"), "shap"), 2), "\n")
for (m in unique(imp_node$model)) cat(m, "proteins:", round(stab(imp_node |> filter(model == m, type == "gene"), "importance"), 2),
                                      " pathways:", round(stab(imp_node |> filter(model == m, type != "gene"), "importance"), 2), "\n")

# top-10 proteins by each method (rep 1, averaged over folds), with Reactome-free sanity check: univariate t
suppressMessages(library(binnr)); aki <- hartman_aki()
tt <- apply(aki$proteomics, 2, function(v) abs(t.test(v ~ aki$design$group)$statistic))
top <- bind_rows(
  shap_node |> filter(source_layer == 0) |> group_by(base) |> summarise(score = mean(shap)) |> mutate(method = "binn SHAP"),
  imp_node |> filter(type == "gene") |> group_by(model, base) |> summarise(score = mean(importance), .groups = "drop") |> rename(method = model)
) |> group_by(method) |> slice_max(score, n = 10) |> mutate(rank = row_number(), univariate_rank = rank(-tt)[base]) |> ungroup()
print(top |> select(method, rank, protein = base, univariate_rank) |> pivot_wider(names_from = method, values_from = c(protein, univariate_rank)), width = 200, n = 10)
saveRDS(list(perf = perf, agree = agree, top = top), file.path(out, "explain_compare.rds"))
