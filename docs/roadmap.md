# binnr: design notes and roadmap

Working notes for contributors. The package is deliberately a *thin* layer:
wherever an established R package already solves a sub-problem well, binnr
should call it rather than re-implement it. This file records what is
delegated, what is still hand-rolled, and why.

## What binnr delegates (as of 0.0.0.9000)

| Concern | Delegated to | Notes |
|---|---|---|
| Graph storage, DAG checks, topological order, reachability, pruning | **igraph** | A `bn_graph` *is* an igraph object with node attributes `type`, `level`, `label`; `bn_nodes()`/`bn_edges()` give tibbles, `as.igraph()`/tidygraph work directly. |
| Degree-preserving and G(n, m) null graphs | **igraph** (`rewire(keeping_degseq())`, `sample_bipartite_gnm()`) | `bn_rewire()` only restricts these to the gene-to-pathway subgraph and keeps the hierarchy. |
| Formula / covariate plumbing | **hardhat** (`mold()`, `forge()`) | Factor levels, missing columns and novel levels are checked by hardhat at prediction time; binnr only centres and scales. |
| Baseline hazard for survival curves | **survival** (`coxph(~ offset(lp))` + `basehaz(centered = FALSE)`) | Breslow ties, matching the Cox partial likelihood used in training. |
| Evaluation metrics (benchmark script and vignette) | **yardstick**, `survival::concordance` | yardstick's `brier_class()` is half the multiclass Brier score. |
| Quasi-variance intervals for method comparisons | **qvcalc** | See `vignettes/benchmark.Rmd`. |
| Tensors, autograd, optimisers | **torch** for R | |

## What is still hand-rolled, and why

- **The message-passing schedule** (`compile_network()`, `bn_module`). Nodes
  are updated once, in topological order, with incoming edges grouped by the
  level of their source so that messages can be gathered with `index_add`
  without concatenating node states. No R GNN library offers a DAG schedule
  with per-edge weights and per-level weight matrices; `torchgnn` provides
  GCN/GraphSAGE-style layers that update *all* nodes at every step, which is a
  different (and for a pathway hierarchy, wasteful) computation. If `torchgnn`
  or a successor gains a scheduled/heterogeneous message-passing layer, the
  `propagate()` method of `bn_module` is the only place that would change.
- **The Cox partial likelihood loss** (`cox_partial_loss()`), written in torch
  ops (log-cumsum-exp over risk sets, Breslow ties). Tested against `coxph`.
- **The training loop** (`fit_bn_impl()`): AdamW, early stopping on a
  validation split, L1 penalty, per-epoch `gc()` to keep R torch memory in
  check.

## Roadmap

1. **luz for the training loop.** `luz` (the Keras-like fitting API for torch
   for R) would replace `fit_bn_impl()`'s loop with `luz::fit()`, giving
   callbacks (early stopping, LR schedules, logging, checkpoints) and metrics
   for free. Blockers: the Cox loss needs the whole batch's risk set (so
   batching must be full-batch or stratified), the L1 penalty and the
   per-epoch `gc()` need custom callbacks, and `bn_tune()` relies on the
   recorded best validation loss. Plan: keep `bn_module` as a plain
   `nn_module`, add a `luz_module()` wrapper and a `bn_control(engine =
   "luz")` option, then switch the default once behaviour is identical in the
   tests.
2. **parsnip engine.** A `binn()` model spec with `set_engine("binnr")` so
   that binnr fits drop into tidymodels workflows and `tune`. hardhat is
   already used, so `mold()`/`forge()` blueprints can be reused.
3. **Cached data download instead of bundled data.** `brca`/`prad` push the
   installed size above 5 MB. Replace with a `tcga_cohort("BRCA")` downloader
   (curatedTCGAData via ExperimentHub, cached with `tools::R_user_dir()`),
   keeping only a small toy cohort in `data/`.
4. **Missing modalities.** Allow a sample to lack an omics layer entirely
   (mask its input nodes) rather than requiring complete cases.
5. **Learned or soft edges.** Treat the knowledge graph as a prior on the
   mask rather than a hard constraint (e.g. sparse learned adjacency with a
   penalty towards the prior), so that "how much of the graph is used" can be
   read off the fit.
6. **Bioconductor submission.** `MultiAssayExperiment` input is supported;
   still needed: BiocCheck compliance, `BiocStyle` vignettes, and ExperimentHub
   data (see item 3).
