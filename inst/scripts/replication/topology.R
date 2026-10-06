# How different are "random" knowledge graphs from the real one?
# Compares binnr's null models (degree = Maslov-Sneppen swaps; pathway_size =
# Caranzano et al.'s column shuffle; bernoulli = their global G(n,m); genes =
# label permutation) on three real graphs, by the quantities that bound what a
# BINN can learn from a gene: direct degree, number of root-to-gene paths,
# reachable roots, and co-membership among genes.
suppressPackageStartupMessages({
  library(binnr)
  library(dplyr); library(tidyr); library(igraph)
})
set.seed(1)
paths_to_roots <- function(g) {
  # number of distinct directed paths from each gene to any root pathway
  ord <- as.integer(topo_sort(g, mode = "in"))   # roots first
  n <- vcount(g); np <- numeric(n)
  roots <- which(degree(g, mode = "out") == 0)
  np[roots] <- 1
  adj <- adjacent_vertices(g, ord, mode = "out")
  for (i in seq_along(ord)) {
    v <- ord[i]; if (v %in% roots) next
    np[v] <- sum(np[as.integer(adj[[i]])])
  }
  setNames(np, V(g)$name)
}
reachable_roots <- function(g) {
  roots <- V(g)$name[degree(g, mode = "out") == 0]
  d <- distances(g, v = V(g)[V(g)$type == "gene"], to = roots, mode = "out")
  rowSums(is.finite(d))
}
gene_cojaccard <- function(g, genes) {
  # Jaccard similarity of pathway membership among a gene set (direct annotations)
  e <- bn_edges(g); e <- e[e$from %in% genes, ]
  sets <- split(e$to, e$from)
  if (length(sets) < 2) return(NA_real_)
  pairs <- combn(names(sets), 2)
  mean(apply(pairs, 2, function(p) length(intersect(sets[[p[1]]], sets[[p[2]]])) / length(union(sets[[p[1]]], sets[[p[2]]]))))
}
summarise_graph <- function(g, label, top_genes) {
  n <- bn_nodes(g); genes <- n$name[n$type == "gene"]
  deg <- degree(g, v = genes, mode = "out")
  np <- paths_to_roots(g)[genes]
  rr <- reachable_roots(g)
  comp <- components(g, mode = "weak")$no
  tg <- intersect(top_genes, genes)
  tibble(
    graph = label, genes = length(genes), pathways = sum(n$type == "pathway"),
    edges = ecount(g), depth = max(n$level), roots = sum(degree(g, mode = "out") == 0),
    gene_deg_median = median(deg), gene_deg_max = max(deg), gene_deg_gini = ineq(deg),
    paths_median = median(np), paths_max = max(np), genes_no_path = sum(np == 0),
    roots_reached_median = median(rr), roots_reached_max = max(rr),
    components = comp,
    top_genes_deg = mean(deg[tg]), top_genes_paths = median(np[tg]),
    top_genes_cojaccard = gene_cojaccard(g, tg)
  )
}
ineq <- function(x) { x <- sort(x); n <- length(x); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) }

run <- function(g, name, top_genes, seeds = 1:5) {
  real <- summarise_graph(g, "Reactome", top_genes)
  nulls <- bind_rows(lapply(c("degree", "pathway_size", "bernoulli", "genes"), function(m) {
    bind_rows(lapply(seeds, function(s) summarise_graph(bn_rewire(g, m, seed = s), m, top_genes) |> mutate(seed = s)))
  }))
  # correlation of per-gene path counts with the real graph (seed 1)
  np_real <- paths_to_roots(g)
  cors <- sapply(c("degree", "pathway_size", "bernoulli", "genes"), function(m) {
    gn <- bn_rewire(g, m, seed = 1); np <- paths_to_roots(gn)
    common <- intersect(names(np_real)[bn_nodes(g)$type == "gene"], names(np))
    cor(np_real[common], np[common], method = "spearman")
  })
  list(dataset = name, real = real, nulls = nulls, path_cor = cors)
}

# 1. binnr's BRCA graph (1,000 genes); "top genes" = the 22 drivers
brca <- tcga_cohort("BRCA")
g1 <- bn_prune(pathway_graph(reactome_tcga$gene_sets, reactome_tcga$hierarchy), genes = colnames(brca$omics$rna))
drivers <- c("TP53", "PIK3CA", "CDH1", "GATA3", "MAP3K1", "KMT2C", "PTEN", "AKT1", "ERBB2", "ESR1", "BRCA1", "BRCA2", "RB1", "NF1", "CBFB", "RUNX1", "TBX3", "FOXA1", "CCND1", "MYC", "MAP2K4", "CDKN1B")
r1 <- run(g1, "TCGA BRCA (binnr, 1,000 genes)", drivers)
# 2. P-NET's 9,229 genes on current Reactome; top genes = P-NET's reported top genes
pc <- pnet_prostate()
g2 <- bn_prune(reactome_graph(), genes = colnames(pc$omics$mut))
pnet_top <- c("AR", "TP53", "PTEN", "RB1", "MDM4", "FGFR1", "MAML3", "PDGFA", "NOTCH1", "EIF3E")
r2 <- run(g2, "P-NET prostate (9,229 genes)", pnet_top)
# 3. Hartman proteomics (UniProt)
aki <- hartman_aki()
g3 <- bn_prune(reactome_graph(identifiers = "uniprot"), genes = colnames(aki$proteomics))
# top proteins: largest |t| between groups
tt <- apply(aki$proteomics, 2, function(v) t.test(v ~ aki$design$group)$statistic)
r3 <- run(g3, "Hartman AKI proteomics (461 proteins)", names(sort(abs(tt), decreasing = TRUE))[1:20])
saveRDS(list(r1, r2, r3), "Sys.getenv("BINNR_REPLICATION_OUT", "replication-out")/topology.rds")
for (r in list(r1, r2, r3)) {
  cat("\n==", r$dataset, "\n")
  print(bind_rows(r$real, r$nulls |> group_by(graph) |> summarise(across(-seed, mean), .groups = "drop")) |>
          select(graph, edges, gene_deg_median, gene_deg_max, gene_deg_gini, paths_median, paths_max, genes_no_path,
                 roots_reached_median, components, top_genes_deg, top_genes_paths, top_genes_cojaccard), width = 200)
  cat("Spearman(per-gene path count, real):", paste(names(r$path_cor), round(r$path_cor, 2), collapse = "  "), "\n")
}
