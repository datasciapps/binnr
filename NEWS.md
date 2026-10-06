# binnr 0.1.0

* The TCGA example cohorts are no longer shipped with the package. `tcga_cohort("BRCA")`
  / `tcga_cohort("PRAD")` build them from `curatedTCGAData` on first use (the same
  derivation as before, now in `R/tcga.R`) and cache the result in
  `tcga_cache_dir()` (`tools::R_user_dir("binnr", "cache")`, or `BINNR_CACHE`);
  `tcga_cached()` tests for a cached copy. Examples, tests and vignettes that need
  the data are skipped when it is not cached. `reactome_tcga` stays bundled.
  The installed package shrinks from about 6 MB to under 1 MB.
* New dataset helpers, all downloaded once into `binnr_cache_dir()`:
  `hartman_aki()` (septic-AKI plasma proteomics shipped with the Python `binn`
  package, 197 x 554 UniProt accessions) and `pnet_prostate()` (P-NET's 1,013
  prostate tumours with mutation and copy number for 9,229 genes, from the
  authors' Zenodo archive, including their train/validation/test split).
* `reactome_graph(identifiers = "uniprot")` builds the knowledge graph from
  `UniProt2Reactome.txt`, for proteomics inputs.
* Continuous integration: R CMD check on Linux/macOS/Windows (with the TCGA
  cohorts built and cached on Linux so vignettes and examples run), test
  coverage (codecov) and a pkgdown site. 176 testthat expectations.

# binnr 0.0.0.9000

* First alpha release for AI4Sci 2026.
* A `bn_graph` is an igraph object (node attributes `type`, `level`,
  `label`), so igraph and tidygraph functions apply directly; null models use
  `igraph::rewire()` and `igraph::sample_bipartite_gnm()`.
* Covariates are handled with hardhat (`mold()`/`forge()`); the baseline
  hazard comes from `survival::basehaz()`.
* Knowledge graphs: `bn_graph()`, `pathway_graph()`, `reactome_graph()`,
  `bn_prune()`, `bn_rewire()` (degree-preserving, pathway-size, Bernoulli and
  gene-permutation null models), `bn_drop_edges()` (incomplete
  knowledge graphs), `bn_pad()` (layered BINN graphs).
* Models: `bn_fit()` with formula, matrix/list and `MultiAssayExperiment`
  interfaces; classification, regression and Cox survival; multi-omics and
  clinical covariates; `bn_arch()` presets `"binn"` and `"gnn"`.
* Methods: `predict()`, `tidy()`, `glance()`, `autoplot()`,
  `node_importance()`, `node_activations()`, `bn_importance_stability()`;
  `update()` refits with a new graph or architecture.
* Example data `brca` (TCGA BRCA) and `reactome_tcga`.
