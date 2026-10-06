#' Reactome pathways for the example data
#'
#' Gene-set membership, pathway hierarchy and pathway names from Reactome,
#' restricted to the genes of the two TCGA example cohorts (see
#' [tcga_cohort()]) and their ancestor pathways.
#' Use [pathway_graph()] to turn it into a knowledge graph, or
#' [reactome_graph()] for the full, current Reactome release.
#'
#' @format A list of three tibbles:
#' \describe{
#'   \item{gene_sets}{`pathway` (Reactome stable id) and `gene` (symbol);
#'     membership at all levels, as in `ReactomePathways.gmt`.}
#'   \item{hierarchy}{`parent` and `child` pathway ids.}
#'   \item{pathways}{`pathway` id and human-readable `name`.}
#' }
#' The Reactome release is stored in `attr(reactome_tcga, "version")`.
#' @source <https://reactome.org/download-data> (CC BY 4.0).
#' @examples
#' g <- pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy,
#'                    labels = reactome_tcga$pathways)
#' g
"reactome_tcga"
