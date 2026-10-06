#' Local cache for downloaded datasets
#'
#' Datasets used in the documentation and benchmarks are not shipped with
#' binnr; they are downloaded on first use and kept here. The directory is
#' the `BINNR_CACHE` environment variable if set, otherwise
#' `tools::R_user_dir("binnr", "cache")`.
#'
#' @return A path.
#' @seealso [tcga_cohort()], [hartman_aki()], [pnet_prostate()],
#'   [reactome_graph()]
#' @export
binnr_cache_dir <- function() {
  env <- Sys.getenv("BINNR_CACHE", unset = "")
  if (nzchar(env)) env else tools::R_user_dir("binnr", "cache")
}

cached_download <- function(url, file, cache_dir = binnr_cache_dir(), refresh = FALSE,
                            what = basename(url)) {
  path <- file.path(cache_dir, file)
  if (refresh || !file.exists(path)) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    cli::cli_inform("Downloading {.file {what}} to {.path {cache_dir}}.")
    withr::local_options(timeout = max(1800, getOption("timeout")))
    utils::download.file(url, path, quiet = TRUE, mode = "wb")
  }
  path
}

#' Septic AKI plasma proteomics (Hartman et al., 2023)
#'
#' The plasma proteomics matrix that ships with the Python `binn` package
#' (Hartman et al., 2023, *Nature Communications*): 554 proteins (UniProt
#' accessions) quantified by data-independent acquisition mass spectrometry
#' in 197 patients with sepsis, labelled by acute kidney injury group. The
#' two files (`sample_datamatrix.csv`, `sample_design_matrix.tsv`) are
#' downloaded from the package's GitHub repository and cached.
#'
#' About 37% of intensities are missing (not quantified). `binn`'s own data
#' loader replaces them with 0 before standardising; `impute = "zero"`
#' reproduces that, `impute = "min"` uses half the protein's minimum observed
#' intensity (a common proteomics convention) and `impute = "none"` returns
#' `NA`s for you to handle. Intensities are log2 scale.
#'
#' The design file labels the groups `1` (n = 74) and `2` (n = 123) without
#' further annotation; see the paper for the clinical definition.
#'
#' Use with [reactome_graph()] and `identifiers = "uniprot"` to obtain a
#' knowledge graph keyed by the same accessions.
#'
#' @param impute How to treat missing intensities: `"zero"`, `"min"` or
#'   `"none"`.
#' @param cache_dir,refresh As for [tcga_cohort()].
#' @return A list with `proteomics` (197 x 554 matrix, samples in rows, UniProt
#'   accessions in columns) and `design` (a tibble with `sample` and `group`,
#'   a factor with levels `"1"` and `"2"`).
#' @references Hartman, E. et al. (2023). Interpreting biologically informed
#'   neural networks for enhanced proteomic biomarker discovery and pathway
#'   analysis. *Nature Communications* 14, 5359.
#'   <https://github.com/InfectionMedicineProteomics/BINN>
#' @examplesIf interactive()
#' aki <- hartman_aki()
#' g <- reactome_graph(identifiers = "uniprot")
#' fit <- bn_fit(aki$proteomics, aki$design$group, graph = g)
#' @export
hartman_aki <- function(impute = c("zero", "min", "none"), cache_dir = binnr_cache_dir(),
                        refresh = FALSE) {
  impute <- rlang::arg_match(impute)
  base <- "https://raw.githubusercontent.com/InfectionMedicineProteomics/BINN/main/binn/data/"
  dm <- cached_download(paste0(base, "sample_datamatrix.csv"), "hartman_aki/sample_datamatrix.csv",
                        cache_dir, refresh)
  de <- cached_download(paste0(base, "sample_design_matrix.tsv"), "hartman_aki/sample_design_matrix.tsv",
                        cache_dir, refresh)
  x <- utils::read.csv(dm, check.names = FALSE)
  design <- utils::read.delim(de)
  mat <- t(as.matrix(x[, -1]))
  colnames(mat) <- x$Protein
  mat <- mat[match(design$sample, rownames(mat)), , drop = FALSE]
  if (impute == "zero") {
    mat[is.na(mat)] <- 0
  } else if (impute == "min") {
    mins <- apply(mat, 2, min, na.rm = TRUE) / 2
    idx <- which(is.na(mat), arr.ind = TRUE)
    mat[idx] <- mins[idx[, 2]]
  }
  list(
    proteomics = mat,
    design = tibble::tibble(sample = design$sample, group = factor(design$group))
  )
}

#' P-NET prostate cancer cohort (Elmarakeby et al., 2021)
#'
#' The processed inputs of P-NET (Elmarakeby et al., 2021, *Nature*): 1,013
#' prostate tumours from Armenia et al. (2018), with somatic mutations and
#' copy number, labelled primary (0) or metastatic (1). The data are the
#' authors' own release on Zenodo (record 10774954, `_database.zip`, 356 MB,
#' AGPL-3.0), downloaded once and cached; only the files needed are read.
#'
#' P-NET codes each gene by three binary inputs: mutated, deep deletion
#' (GISTIC -2) and high amplification (GISTIC 2). `cna = "binary"` reproduces
#' that as two layers `cna_del` and `cna_amp`; `cna = "gistic"` keeps the
#' thresholded -2..2 values in one layer. Genes follow P-NET's configuration
#' (`onsplit_average_reg_10_tanh_large_testing.py`): the union of genes
#' present in either matrix, restricted to P-NET's list of genes expressed in
#' TCGA prostate cancer or known cancer genes and to HGNC protein-coding
#' genes, with genes absent from a matrix filled with 0. This gives the 9,229
#' genes of the paper. `genes = "all"` keeps every gene in the matrices and a
#' character vector restricts to those symbols.
#'
#' @param cna How to encode copy number: `"binary"` (P-NET) or `"gistic"`.
#' @param genes `"pnet"` for the archive's gene list, `"all"` for every gene in
#'   the matrices, or a character vector of symbols.
#' @param cache_dir,refresh As for [tcga_cohort()].
#' @return A list with `omics` (named list of samples x genes 0/1 or integer
#'   matrices), `clinical` (tibble: `sample`, `response` factor
#'   `"primary"`/`"metastatic"`) and `splits` (tibble: `sample`, `set` with
#'   P-NET's `train`/`validation`/`test` assignment, where available).
#' @references Elmarakeby, H. A. et al. (2021). Biologically informed deep
#'   neural network for prostate cancer discovery. *Nature* 598, 348-352.
#'   Data: \doi{10.5281/zenodo.10774954}.
#' @examplesIf interactive()
#' pc <- pnet_prostate()
#' str(pc$omics, max.level = 1)
#' table(pc$clinical$response, pc$splits$set)
#' @export
pnet_prostate <- function(cna = c("binary", "gistic"), genes = "pnet",
                          cache_dir = binnr_cache_dir(), refresh = FALSE) {
  cna <- rlang::arg_match(cna)
  dir <- file.path(cache_dir, "pnet")
  files <- c(
    mut = "_database/prostate/processed/P1000_final_analysis_set_cross_important_only.csv",
    cna = "_database/prostate/processed/P1000_data_CNA_paper.csv",
    response = "_database/prostate/processed/response_paper.csv",
    genes = "_database/genes/tcga_prostate_expressed_genes_and_cancer_genes.csv",
    coding = "_database/genes/HUGO_genes/protein-coding_gene_with_coordinate_minimal.txt",
    train = "_database/prostate/splits/training_set_0.csv",
    validation = "_database/prostate/splits/validation_set.csv",
    test = "_database/prostate/splits/test_set.csv"
  )
  if (refresh || !all(file.exists(file.path(dir, files)))) {
    zip <- cached_download("https://zenodo.org/records/10774954/files/_database.zip?download=1",
                           "pnet/_database.zip", cache_dir, refresh, what = "_database.zip (356 MB)")
    utils::unzip(zip, files = unname(files), exdir = dir, overwrite = TRUE)
  }
  read <- function(key, ...) utils::read.csv(file.path(dir, files[[key]]), check.names = FALSE, ...)
  mut <- read("mut", row.names = 1)
  cn <- read("cna", row.names = 1)
  response <- read("response")
  samples <- intersect(intersect(rownames(mut), rownames(cn)), response$id)
  gene_univ <- union(colnames(mut), colnames(cn))
  if (identical(genes, "pnet")) {
    coding <- utils::read.delim(file.path(dir, files[["coding"]]), header = FALSE)[[4]]
    genes <- intersect(intersect(gene_univ, read("genes")$genes), coding)
  } else if (identical(genes, "all")) {
    genes <- gene_univ
  } else {
    genes <- intersect(gene_univ, as.character(genes))
  }
  genes <- sort(genes)
  # genes missing from a matrix are absent because no sample carries the event
  fill <- function(df) {
    out <- matrix(0L, length(samples), length(genes), dimnames = list(samples, genes))
    hit <- intersect(genes, colnames(df))
    out[, hit] <- as.integer(as.matrix(df[samples, hit]))
    out
  }
  m <- fill(mut); m <- (m > 0) * 1L
  cnm <- fill(cn)
  omics <- if (cna == "binary") {
    list(mut = m, cna_del = (cnm <= -2) * 1L, cna_amp = (cnm >= 2) * 1L)
  } else {
    list(mut = m, cna = cnm)
  }
  omics <- lapply(omics, function(x) { dimnames(x) <- list(samples, genes); x })
  resp <- response$response[match(samples, response$id)]
  splits <- do.call(rbind, lapply(c("train", "validation", "test"), function(s) {
    ids <- read(s)$id
    data.frame(sample = ids, set = s)
  }))
  structure(
    list(
      omics = omics,
      clinical = tibble::tibble(sample = samples,
                                response = factor(resp, 0:1, c("primary", "metastatic"))),
      splits = tibble::as_tibble(splits[splits$sample %in% samples, ])
    ),
    source = "P-NET _database.zip, Zenodo 10.5281/zenodo.10774954 (Armenia et al. 2018 cohort)"
  )
}
