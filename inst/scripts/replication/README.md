# Replication scripts

Anchoring binnr against prior work (see `docs/` and the 2 Oct 2026 report):

- `run_binn.py` — Hartman et al. (2023) septic-AKI proteomics with the Python
  `binn` package (v0.1.1, its defaults), 5x5 CV; exports folds, predictions,
  the layered connectivity graph and SHAP importances.
- `compare.R` — binnr on the identical folds: binn's exact layered graph
  (stratified) and binnr's topological Reactome DAG, plus a lasso.
- `explain_compare.R` — agreement and stability of explanations (binn SHAP vs
  binnr gradient x activation).
- `topology.R` — what binnr's null models do to gene degree, path counts and
  co-membership on three real graphs (BRCA, P-NET 9,229 genes, Hartman).
- `nulls.R` — real vs null graphs on Hartman AKI (3x5 CV) and on P-NET's own
  split (5 seeds).

Outputs go to `$BINNR_REPLICATION_OUT` (default `replication-out/`). Data come
from `hartman_aki()`, `pnet_prostate()`, `tcga_cohort()`.
