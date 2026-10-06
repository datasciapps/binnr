# Per-gene path budgets on P-NET's gene universe (current Reactome) under the
# real graph and two null models; informative genes chosen from the data
# (P-NET training split, Fisher test of any alteration vs metastasis).
suppressMessages({ devtools::load_all(".", quiet = TRUE); library(igraph) })
paths_to_roots <- function(g) {
  ord <- as.integer(topo_sort(g, mode = "in")); n <- vcount(g); np <- numeric(n)
  roots <- which(degree(g, mode = "out") == 0); np[roots] <- 1
  adj <- adjacent_vertices(g, ord, mode = "out")
  for (i in seq_along(ord)) { v <- ord[i]; if (!(v %in% roots)) np[v] <- sum(np[as.integer(adj[[i]])]) }
  setNames(np, V(g)$name)
}
pc <- pnet_prostate()
g <- bn_prune(reactome_graph(), genes = colnames(pc$omics$mut))
genes <- bn_nodes(g)$name[bn_nodes(g)$type == "gene"]
tr <- pc$splits$sample[pc$splits$set == "train"]
y <- pc$clinical$response[match(tr, pc$clinical$sample)]
alt <- (pc$omics$mut[tr, genes] + pc$omics$cna_del[tr, genes] + pc$omics$cna_amp[tr, genes]) > 0
p <- apply(alt, 2, function(a) if (length(unique(a)) < 2) 1 else fisher.test(table(a, y))$p.value)
top <- names(sort(p))[1:20]
real <- paths_to_roots(g)[genes]
deg <- degree(g, v = genes, mode = "out")
out <- lapply(c("degree", "pathway_size", "genes"), function(m) {
  gn <- bn_rewire(g, m, seed = 1); pn <- paths_to_roots(gn)
  # genes left without any pathway are dropped from the null graph: 0 paths
  dn <- setNames(numeric(length(genes)), genes); kept <- intersect(genes, V(gn)$name)
  dn[kept] <- degree(gn, v = kept, mode = "out")
  pnull <- setNames(numeric(length(genes)), genes); pnull[kept] <- pn[kept]
  data.frame(null = m, gene = genes, real = real, null_paths = pnull,
             deg_real = deg, deg_null = dn, top = genes %in% top)
})
d <- do.call(rbind, out)
# exchangeable genes: identical annotation sets
e <- bn_edges(g); key <- tapply(e$to[e$from %in% genes], e$from[e$from %in% genes], function(s) paste(sort(s), collapse = "|"))
tab <- table(key); exch <- sum(tab[tab > 1]) / length(genes)
saveRDS(list(paths = d, top = top, p = p[top], exch = exch, n_genes = length(genes), g_summary = summary(g)),
        "talk/data/null_paths.rds")
cat("top genes:", paste(top, collapse = ", "), "\nexchangeable:", round(exch, 3), "\n")
print(aggregate(cbind(real, null_paths, deg_real, deg_null) ~ null + top, d, median))
print(aggregate(cbind(real, null_paths, deg_real, deg_null) ~ null, d, max))
