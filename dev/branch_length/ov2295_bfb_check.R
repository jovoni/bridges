# OV2295 (500 kb, wgd = TRUE): BFB-like nodes and both BFB tests on chr 6
# and chr 22, at jitter 0 and 2. Run before and after a change to
# reconstruct_tree() to compare.
#   Rscript dev/branch_length/ov2295_bfb_check.R [output.rds]
suppressMessages(library(bridges))
source("dev/wgd/real_data_prep.R")
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "ov2295_bfb_check.rds"
chrs <- as.character(1:22)
rd <- load_ov2295(1)
rows <- list()
for (jit in c(0, 2)) {
  t0 <- Sys.time()
  res <- fit(rd$data, chromosomes = chrs, alleles = "CN", k_jitter_fix = jit, wgd = TRUE, verbose = FALSE)
  t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  br <- detect_bfb_branches(res)
  cl <- detect_bfb_clades(res, seed = 1)
  n_bfb <- sapply(chrs, function(ch) sum(unlist(res$reconstructions[[ch]][["CN"]]$deltas) > 0, na.rm = TRUE))
  depth <- ape::node.depth.edgelength(res$tree)[seq_along(res$tree$tip.label)]
  wgd_len <- res$tree$edge.length[res$tree$edge[, 2] %in% res$wgd$edges$child]
  rows[[length(rows) + 1]] <- data.frame(
    jitter = jit, chr = chrs, bfb_nodes = n_bfb,
    branches_adj = br$adj.pval[match(chrs, br$chr)], clades_adj = cl$adj.pval[match(chrs, cl$chr)],
    clade_cells = cl$clade_n_cells[match(chrs, cl$chr)],
    fit_sec = t_fit, mean_depth = mean(depth), max_wgd_edge = max(wgd_len))
  cat(sprintf("jitter %d: fit %.1f s, mean root-to-tip %.1f, WGD edge lengths %s\n",
              jit, t_fit, mean(depth), paste(round(sort(wgd_len, decreasing = TRUE)[1:3], 1), collapse = ", ")))
}
d <- do.call(rbind, rows)
saveRDS(d, out_file)
cat("\nchromosomes significant in either test (adj < 0.05):\n")
print(d[d$branches_adj < 0.05 | d$clades_adj < 0.05, c("jitter", "chr", "bfb_nodes", "branches_adj", "clades_adj", "clade_cells")],
      row.names = FALSE, digits = 2)
cat("\nchr 6 and 22:\n")
print(d[d$chr %in% c("6", "22"), c("jitter", "chr", "bfb_nodes", "branches_adj", "clades_adj", "clade_cells")], row.names = FALSE, digits = 2)
