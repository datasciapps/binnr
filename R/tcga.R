#' TCGA example cohorts, downloaded and cached
#'
#' Builds the breast (`"BRCA"`) or prostate (`"PRAD"`) cancer example cohort
#' used throughout the binnr documentation from The Cancer Genome Atlas, via
#' the Bioconductor package `curatedTCGAData` (Firehose run 2016-01-28,
#' `version = "2.1.1"`). The first call downloads the assays (a few hundred
#' megabytes through ExperimentHub) and keeps a compact processed copy in a
#' local cache; later calls read the cache, so the data are not shipped with
#' the package.
#'
#' Each cohort keeps `n_genes` Reactome-annotated genes: a short list of
#' well-known driver genes for that cancer, so that the mutation and copy
#' number layers are informative, plus the most variable genes by RNA-seq.
#' Only primary tumours observed in every omics layer are kept.
#'
#' **BRCA** (960 tumours): `rna` (log2 normalised RSEM), `cnv` (GISTIC
#' thresholded copy number, -2 to 2) and `mut` (1 if the gene carries a
#' non-silent somatic mutation); clinical columns `patient`, `age`, `stage`
#' (AJCC I-IV), `pam50` (Luminal A/B, HER2-enriched, Basal-like; `NA` for
#' Normal-like or missing), `os_time` (days) and `os_event` (1 = died).
#'
#' **PRAD** (491 tumours): `rna` and `cnv` only, since the mutation layer
#' covers a third of the cohort; clinical columns `patient`, `age`, `psa`
#' (pre-operative, ng/mL), `t_stage`, `gleason` (6-10), `grade` (`"high"` for
#' Gleason 8-10, as in P-NET) and `recurrence`.
#'
#' @param cohort `"BRCA"` or `"PRAD"`.
#' @param cache_dir Directory for the processed cohorts; see [binnr_cache_dir()].
#' @param refresh Rebuild the cohort even if a cached copy exists.
#' @param n_genes Number of genes to keep.
#' @param gene_sets A data frame of Reactome `pathway`/`gene` membership used
#'   to restrict the gene universe. The default, [reactome_tcga], reproduces the
#'   published example data; pass the `gene_sets` element of
#'   [reactome_graph()]'s source files to select from all of Reactome.
#' @param quiet Suppress progress messages.
#' @return A list with elements `omics` (named list of patients x genes
#'   matrices) and `clinical` (a tibble with one row per patient), with
#'   attribute `"source"`.
#' @seealso [tcga_cached()] to test for a cached copy without downloading;
#'   [reactome_tcga]; [mae_to_omics()] to use any `MultiAssayExperiment`
#'   directly.
#' @examplesIf binnr::tcga_cached("BRCA")
#' brca <- tcga_cohort("BRCA")
#' str(brca$omics, max.level = 1)
#' table(brca$clinical$pam50, useNA = "ifany")
#' @export
tcga_cohort <- function(cohort = c("BRCA", "PRAD"), cache_dir = tcga_cache_dir(),
                        refresh = FALSE, n_genes = 1000L,
                        gene_sets = binnr::reactome_tcga$gene_sets, quiet = FALSE) {
  cohort <- rlang::arg_match(cohort)
  path <- tcga_cache_path(cohort, cache_dir)
  if (!refresh && file.exists(path)) {
    return(readRDS(path))
  }
  rlang::check_installed(
    c("curatedTCGAData", "MultiAssayExperiment", "SummarizedExperiment",
      "TCGAutils", "S4Vectors"),
    reason = "to download TCGA cohorts (Bioconductor)."
  )
  if (!quiet) {
    cli::cli_inform(c(
      "Downloading TCGA {cohort} via {.pkg curatedTCGAData}; this happens once.",
      i = "The processed cohort is cached in {.path {cache_dir}}."
    ))
  }
  out <- switch(cohort,
    BRCA = build_brca(n_genes, gene_sets),
    PRAD = build_prad(n_genes, gene_sets)
  )
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, path, compress = "xz")
  out
}

#' @rdname tcga_cohort
#' @return `tcga_cached()` returns `TRUE` if a processed copy of the cohort
#'   is in the cache (so [tcga_cohort()] will not download anything).
#' @export
tcga_cached <- function(cohort = c("BRCA", "PRAD"), cache_dir = tcga_cache_dir()) {
  cohort <- rlang::arg_match(cohort)
  file.exists(tcga_cache_path(cohort, cache_dir))
}

#' @rdname tcga_cohort
#' @export
tcga_cache_dir <- function() binnr_cache_dir()

tcga_cache_path <- function(cohort, cache_dir) {
  file.path(cache_dir, paste0("tcga_", tolower(cohort), ".rds"))
}

# Aligned omics matrices (patients x genes) and colData for one cohort -------
tcga_omics <- function(cohort, layers, drivers, n_genes, gene_sets) {
  patient <- function(barcode) substr(barcode, 1L, 12L)
  assays <- c(rna = "RNASeq2GeneNorm", cnv = "GISTIC_ThresholdedByGene",
              mut = "Mutation")[layers]
  mae <- curatedTCGAData::curatedTCGAData(cohort, assays = unname(assays),
                                          version = "2.1.1", dry.run = FALSE)
  mae <- TCGAutils::TCGAprimaryTumors(mae)
  mae <- MultiAssayExperiment::intersectColumns(mae)
  nm_of <- function(a) grep(paste0("_", a, "-"), names(mae), value = TRUE)

  rna <- SummarizedExperiment::assay(mae[[nm_of(assays[["rna"]])]])
  colnames(rna) <- patient(colnames(rna))
  cnv_se <- mae[[nm_of(assays[["cnv"]])]]
  cnv <- SummarizedExperiment::assay(cnv_se)
  storage.mode(cnv) <- "integer"
  rownames(cnv) <- SummarizedExperiment::rowData(cnv_se)$Gene.Symbol
  cnv <- cnv[!duplicated(rownames(cnv)), ]
  colnames(cnv) <- patient(colnames(cnv))
  samples <- list(colnames(rna), colnames(cnv))

  mut_list <- NULL
  if ("mut" %in% layers) {
    gr <- methods::as(mae[[nm_of(assays[["mut"]])]], "GRangesList")
    mut_list <- lapply(gr, function(g) {
      keep <- !S4Vectors::mcols(g)$Variant_Classification %in% c("Silent", "RNA")
      unique(S4Vectors::mcols(g)$Hugo_Symbol[keep])
    })
    names(mut_list) <- patient(names(mut_list))
    samples <- c(samples, list(names(mut_list)))
  }
  samples <- sort(unique(Reduce(intersect, samples)))

  # Gene selection: Reactome-annotated, fully observed RNA, present in CNV;
  # drivers first, then the most variable genes
  genes <- intersect(intersect(rownames(rna), rownames(cnv)), unique(gene_sets$gene))
  rna <- rna[genes, samples]
  genes <- genes[rowSums(is.na(rna)) == 0]
  v <- apply(rna[genes, ], 1, stats::var)
  drivers <- intersect(drivers, genes)
  top <- setdiff(names(sort(v, decreasing = TRUE)), drivers)
  genes <- sort(c(drivers, top[seq_len(max(0L, n_genes - length(drivers)))]))

  omics <- list(rna = t(round(rna[genes, samples], 2)), cnv = t(cnv[genes, samples]))
  if (!is.null(mut_list)) {
    mut <- t(vapply(mut_list[samples], function(g) genes %in% g, logical(length(genes))))
    dimnames(mut) <- list(samples, genes)
    storage.mode(mut) <- "integer"
    omics$mut <- mut
  }
  cd <- as.data.frame(MultiAssayExperiment::colData(mae))
  cd <- cd[match(samples, cd$patientID), ]
  list(omics = omics, cd = cd, samples = samples, genes = genes)
}

build_brca <- function(n_genes, gene_sets) {
  b <- tcga_omics("BRCA", c("rna", "cnv", "mut"), n_genes = n_genes, gene_sets = gene_sets,
                  drivers = c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C",
                              "PTEN", "AKT1", "ERBB2", "ESR1", "BRCA1", "BRCA2", "RB1",
                              "NF1", "CBFB", "RUNX1", "TBX3", "FOXA1", "CCND1", "MYC",
                              "MAP2K4", "CDKN1B"))
  cd <- b$cd
  stage <- toupper(sub("^stage ([ivx]+).*$", "\\1", cd$pathologic_stage))
  stage[!stage %in% c("I", "II", "III", "IV")] <- NA
  pam50 <- cd$PAM50.mRNA
  pam50[pam50 == "Normal-like"] <- NA # only 8 samples
  time <- ifelse(is.na(cd$days_to_death), cd$days_to_last_followup, cd$days_to_death)
  clinical <- tibble::tibble(
    patient = b$samples,
    age = as.integer(cd$years_to_birth),
    stage = factor(stage, levels = c("I", "II", "III", "IV")),
    pam50 = factor(pam50, levels = c("Luminal A", "Luminal B", "HER2-enriched",
                                     "Basal-like")),
    os_time = as.numeric(time),
    os_event = as.integer(cd$vital_status)
  )
  clinical$os_time[!is.na(clinical$os_time) & clinical$os_time <= 0] <- NA
  structure(list(omics = b$omics, clinical = clinical),
            source = "curatedTCGAData 2.1.1 (TCGA BRCA, Firehose 2016-01-28)")
}

build_prad <- function(n_genes, gene_sets) {
  # RNA and copy number only: the mutation layer covers just 332 of the ~490
  # primary tumours with RNA and CNV. The outcome of interest is Gleason grade
  # (>= 8 vs <= 7), as in P-NET (Elmarakeby et al., 2021).
  p <- tcga_omics("PRAD", c("rna", "cnv"), n_genes = n_genes, gene_sets = gene_sets,
                  drivers = c("TP53", "PTEN", "SPOP", "FOXA1", "AR", "ERG", "TMPRSS2",
                              "RB1", "MYC", "CDK12", "KMT2C", "KMT2D", "ATM", "BRCA2",
                              "APC", "CTNNB1", "PIK3CA", "ZBTB16", "NKX3-1", "IDH1"))
  cd <- p$cd
  gleason <- as.integer(cd$patient.stage_event.gleason_grading.gleason_score)
  psa <- suppressWarnings(as.numeric(cd$patient.clinical_cqcf.psa_result_preop))
  t_stage <- toupper(sub("^(t[1-4]).*$", "\\1", cd$pathology_T_stage))
  t_stage[!t_stage %in% c("T1", "T2", "T3", "T4")] <- NA
  clinical <- tibble::tibble(
    patient = p$samples,
    age = as.integer(cd$years_to_birth),
    psa = psa,
    t_stage = factor(t_stage, levels = c("T1", "T2", "T3", "T4")),
    gleason = gleason,
    grade = factor(ifelse(gleason >= 8L, "high", "low"), levels = c("low", "high")),
    recurrence = factor(cd$patient.biochemical_recurrence, levels = c("no", "yes"))
  )
  structure(list(omics = p$omics, clinical = clinical),
            source = "curatedTCGAData 2.1.1 (TCGA PRAD, Firehose 2016-01-28)")
}
