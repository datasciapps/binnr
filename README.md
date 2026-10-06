# binnr

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/datasciapps/binnr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/datasciapps/binnr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**binnr** fits *biologically informed neural networks* (BINNs) in R. Omics
features feed genes, genes feed Reactome (or any other) pathways, and
pathways feed their parent pathways. It is a thin, idiomatic R front-end
built on plain [torch](https://torch.mlverse.org): the network is a
**message-passing graph neural network** on the knowledge DAG, with no
custom graph library, so:

- the classic feedforward BINN (P-NET, BINN) is a special case, and it is the
  default;
- no dummy padding nodes are needed for ragged hierarchies;
- each BINN assumption (scalar pathway nodes, one weight per edge, readout
  from the top pathways only) is a single argument you can relax;
- multi-omics inputs, clinical covariates, and classification, regression and
  Cox survival outcomes are supported out of the box;
- the graph is data, not architecture: any directed acyclic graph from
  features to concepts can be supplied as an edge list, and null models
  (`bn_rewire()`), incomplete graphs (`bn_drop_edges()`) and padded strata
  (`bn_pad()`) are one-line variations.

It uses standard R and Bioconductor objects (matrices, data frames, formulas,
`survival::Surv`, `MultiAssayExperiment`) and returns tidy output.

> **Status:** alpha / proof of concept. The API may change.

## Installation

```r
# install.packages("pak")
pak::pak("datasciapps/binnr")
torch::install_torch()   # once, to download libtorch
```

## Example

```r
library(binnr)
library(survival)

# Reactome knowledge graph (bundled snapshot; see reactome_graph() for the full release)
g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
                   labels = reactome_tcga$pathways)

# TCGA example cohorts, built from curatedTCGAData (Bioconductor) on first use
# and cached locally; nothing is shipped with the package
brca <- tcga_cohort("BRCA")   # 960 tumours: RNA, copy number, mutations
prad <- tcga_cohort("PRAD")   # 491 tumours: RNA, copy number

# TCGA breast cancer: RNA + copy number + mutations, PAM50 subtype
fit <- bn_fit(pam50 ~ 1, data = brca$clinical, omics = brca$omics, graph = g)
predict(fit, brca$omics, type = "prob")

# TCGA prostate cancer: Gleason grade >= 8 from RNA + copy number (the P-NET task)
fit_prad <- bn_fit(grade ~ 1, data = prad$clinical, omics = prad$omics, graph = g)

# Overall survival with clinical covariates, as a relaxed GNN
fit_os <- bn_fit(Surv(os_time, os_event) ~ age + stage,
                 data = brca$clinical, omics = brca$omics, graph = g,
                 arch = bn_arch("gnn"))

# Which pathways drive the predictions?
node_importance(fit, brca$omics)

# Per-patient pathway activity scores, for standard statistics and plots
act <- node_activations(fit, brca$omics, type = "pathway")

# Is it the biology or just the sparsity? Refit on a degree-preserving null graph,
# or on an incomplete knowledge graph with half of the annotations removed.
fit_null <- update(fit, graph = bn_rewire(g, seed = 1))
fit_half <- update(fit, graph = bn_drop_edges(g, prop = 0.5, seed = 1))
```

## BINNs as GNNs

| BINN assumption | `bn_arch()` argument | relaxation |
|---|---|---|
| one scalar per pathway | `hidden_dim = 1` | vector embeddings |
| a unique weight per edge | `weights = "edge"` | `"shared"` (GCN-style mean aggregation) |
| predict from top pathways | `readout = "roots"` | `"pathways"` |
| stratified schedule, dummy nodes | not needed | topological schedule follows the DAG (`bn_pad()` shows the cost) |

A unit test checks that the default architecture is numerically identical to a
masked dense feedforward BINN.

## Learn more

- `vignette("binnr")`: getting started
- `vignette("benchmark")`: comparison with glmnet, ranger and xgboost on TCGA
  breast and prostate cancer, with corrected resampled *t* intervals and
  paired comparisons against padded, null and incomplete knowledge graphs

## Citation

Selby, D., Pavliuk, D., Grossmann, G. and Vollmer, S. (2026). *Open tools for
more accessible and transparent biologically informed deep learning.*
2nd International Conference on AI for Science (AI4Sci 2026), Mainz.

## Authors

David Selby, Daria Pavliuk, Gerrit Grossmann and Sebastian Vollmer, Deutsches
Forschungszentrum für Künstliche Intelligenz (DFKI).
