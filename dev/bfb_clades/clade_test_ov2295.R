# OV2295 (500 kb, k_jitter_fix = 2): detect_bfb_clades() vs detect_bfb_branches().
# Clone labels are used only to describe the best clades.
#   Rscript dev/bfb_clades/clade_test_ov2295.R
suppressMessages(library(bridges))
source("dev/wgd/real_data_prep.R")
chrs <- as.character(1:22); rd <- load_ov2295(1)
for (w in c(FALSE, TRUE)) {
  res <- fit(rd$data, chromosomes = chrs, alleles = "CN", k_jitter_fix = 2, wgd = w, verbose = FALSE)
  g <- detect_bfb_branches(res)
  t <- system.time(cl <- detect_bfb_clades(res, n_sim = 999, seed = 1))[["elapsed"]]
  cl$global_adj <- g$adj.pval[match(cl$chr, g$chr)]
  cl <- cl[order(cl$p.scan), ]
  cat(sprintf("\n== wgd = %s (clade test %.0f s) ==\n", w, t))
  cat(sprintf("significant: genome-wide test: %s | clade test: %s\n",
              paste(g$chr[g$adj.pval < 0.05], collapse = ",") , paste(cl$chr[cl$adj.pval < 0.05], collapse = ",")))
  top <- head(cl, 5)
  for (i in seq_len(nrow(top))) {
    tb <- table(rd$clones[top$clade_cells[[i]]])
    cat(sprintf("  chr %-2s p.scan %.3g adj %.3g | best clade: %d of %d nodes BFB-like, %d cells (%s) | genome-wide adj %.3g\n",
                top$chr[i], top$p.scan[i], top$adj.pval[i], top$clade_bfb_nodes[i], top$clade_nodes[i],
                top$clade_n_cells[i], paste(names(tb), tb, sep = ":", collapse = " "), top$global_adj[i]))
  }
}
