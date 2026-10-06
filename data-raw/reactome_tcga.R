# Build the bundled `reactome_tcga` snapshot (the TCGA cohorts themselves are
# downloaded on demand by tcga_cohort(); see R/tcga.R).
#
# Source: curatedTCGAData (Bioconductor), Firehose 2016-01-28 run, plus the
# Reactome "current" release downloaded from reactome.org.
#
# To keep the installed package small (< 5 MB), each cohort keeps 1,000
# Reactome-annotated genes: a handful of well-known driver genes for that
# cancer (so that mutation / copy-number layers are informative) plus the most
# variable genes by RNA-seq. Values are RNA-seq (log2 RSEM, as provided),
# GISTIC thresholded copy number (-2..2) and, for BRCA, binary non-silent
# somatic mutation. Only primary tumours observed in every layer are kept.
#
# `reactome_tcga` holds the Reactome gene sets, hierarchy and names restricted
# to the union of genes in both cohorts and their ancestor pathways.

suppressPackageStartupMessages({
  library(curatedTCGAData)
  library(MultiAssayExperiment)
  library(RaggedExperiment)
  library(TCGAutils)
})

n_genes <- 1000L
patient <- function(barcode) substr(barcode, 1L, 12L)

# ---- Reactome -----------------------------------------------------------------
cache <- tempfile("reactome")
dir.create(cache)
base <- "https://reactome.org/download/current/"
utils::download.file(paste0(base, "ReactomePathways.gmt.zip"),
                     file.path(cache, "gmt.zip"), quiet = TRUE)
utils::download.file(paste0(base, "ReactomePathwaysRelation.txt"),
                     file.path(cache, "rel.txt"), quiet = TRUE)
utils::download.file(paste0(base, "ReactomePathways.txt"),
                     file.path(cache, "names.txt"), quiet = TRUE)
reactome_version <- readLines(
  "https://reactome.org/ContentService/data/database/version", warn = FALSE
)

gmt_lines <- readLines(unz(file.path(cache, "gmt.zip"), "ReactomePathways.gmt"))
gmt <- strsplit(gmt_lines, "\t", fixed = TRUE)
gene_sets <- data.frame(
  pathway = rep(vapply(gmt, `[`, "", 2L), lengths(gmt) - 2L),
  gene = unlist(lapply(gmt, `[`, -(1:2)), use.names = FALSE)
)
rel <- utils::read.delim(file.path(cache, "rel.txt"), header = FALSE,
                         col.names = c("parent", "child"))
rel <- rel[startsWith(rel$parent, "R-HSA-") & startsWith(rel$child, "R-HSA-"), ]
nm <- utils::read.delim(file.path(cache, "names.txt"), header = FALSE, quote = "",
                        col.names = c("pathway", "name", "species"))
nm <- nm[nm$species == "Homo sapiens", c("pathway", "name")]
# CRAN/Bioconductor require ASCII data: transliterate (e.g. Greek letters)
nm$name <- stringi::stri_trans_general(nm$name, "Any-Latin; Latin-ASCII")

# ---- TCGA cohorts (not shipped) ----------------------------------------------------
# The BRCA and PRAD example cohorts are built on demand by tcga_cohort() in
# R/tcga.R, from curatedTCGAData, and cached locally. Here they are built
# only to fix the gene universe of the bundled Reactome snapshot.
devtools::load_all(".", quiet = TRUE)
b <- list(genes = colnames(tcga_cohort("BRCA", gene_sets = gene_sets)$omics$rna))
p <- list(genes = colnames(tcga_cohort("PRAD", gene_sets = gene_sets)$omics$rna))

# ---- Reactome restricted to the genes in both cohorts ---------------------------
genes <- union(b$genes, p$genes)
reactome_tcga <- list(
  gene_sets = tibble::as_tibble(gene_sets[gene_sets$gene %in% genes, ]),
  hierarchy = tibble::as_tibble(rel),
  pathways = tibble::as_tibble(nm)
)
# Keep only hierarchy edges between pathways that are ancestors of used pathways
used <- unique(reactome_tcga$gene_sets$pathway)
repeat {
  up <- unique(c(used, rel$parent[rel$child %in% used]))
  if (length(up) == length(used)) break
  used <- up
}
reactome_tcga$hierarchy <- reactome_tcga$hierarchy[
  reactome_tcga$hierarchy$child %in% used & reactome_tcga$hierarchy$parent %in% used, ]
reactome_tcga$pathways <- reactome_tcga$pathways[reactome_tcga$pathways$pathway %in% used, ]
attr(reactome_tcga, "version") <- paste("Reactome", reactome_version)

usethis::use_data(reactome_tcga, overwrite = TRUE, compress = "xz")
