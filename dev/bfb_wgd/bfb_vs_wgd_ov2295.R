# OV2295 (500 kb, k_jitter_fix = 2): BFB detection with and without WGD.
#   Rscript dev/bfb_wgd/bfb_vs_wgd_ov2295.R
suppressMessages(library(bridges))
source("dev/wgd/real_data_prep.R")
chrs <- as.character(1:22)
rd <- load_ov2295(1)

# ---- Tree-based detection ----
fits <- list(no_wgd = fit(rd$data, chromosomes = chrs, alleles = "CN", k_jitter_fix = 2, verbose = FALSE),
             wgd = fit(rd$data, chromosomes = chrs, alleles = "CN", k_jitter_fix = 2, wgd = TRUE, verbose = FALSE))
k_cell <- fits$wgd$wgd$counts
for (nm in names(fits)) {
  res <- fits[[nm]]
  det <- detect_bfb_branches(res)
  det <- det[order(det$adj.pval), ]
  tree <- res$tree; n_tip <- length(tree$tip.label)
  # BFB-like nodes whose two children subtrees mix cells of different WGD counts
  mixed <- 0; total <- 0
  for (chr in chrs) {
    d <- unlist(res$reconstructions[[chr]][["CN"]]$deltas)
    nodes <- as.integer(names(d)[which(d > 0)])
    total <- total + length(nodes)
    for (node in nodes) {
      kids <- tree$edge[tree$edge[, 1] == node, 2]
      ks <- lapply(kids, function(x) {
        tips <- if (x <= n_tip) tree$tip.label[x] else ape::extract.clade(tree, x)$tip.label
        unique(k_cell[tips])
      })
      if (!identical(sort(ks[[1]]), sort(ks[[2]])) || length(unique(unlist(ks))) > 1) mixed <- mixed + 1
    }
  }
  cat(sprintf("\n== %s ==\nBFB-like nodes: %d over all chromosomes, of which %d join cells with different WGD counts\n",
              nm, total, mixed))
  cat(sprintf("chromosomes with adj.pval < 0.05: %d of %d\n", sum(det$adj.pval < 0.05), nrow(det)))
  print(head(data.frame(chr = det$chr, mean = round(det$mean, 3), N_successes = det$N_successes,
                        N_trials = det$N_trials, adj.pval = signif(det$adj.pval, 2)), 8), row.names = FALSE)
}

# ---- Per-cell detection: near-diploid vs doubled cells ----
halved <- rd$data
hi <- halved$cell_id %in% names(k_cell)[k_cell > 0]
halved$CN[hi] <- floor(halved$CN[hi] / 2)
group <- ifelse(k_cell > 0, "doubled", "near-diploid")
for (variant in c("raw", "doubled cells halved")) {
  d <- if (variant == "raw") rd$data else halved
  det <- bfb_detect_batch(d, alleles = "CN", genome = "hg19")
  inf <- det[!is.na(det$admits_bfb), ]
  inf$group <- group[inf$cell_id]
  cat(sprintf("\nper-cell detection (%s): informative count-vectors and fraction admitting BFB\n", variant))
  print(aggregate(admits_bfb ~ group, data = inf, FUN = function(x) c(n = length(x), admits = round(mean(x), 3))))
  det$group <- group[det$cell_id]
  cat("fraction of count-vectors that are informative (not NA):\n")
  print(aggregate(admits_bfb ~ group, data = transform(det, admits_bfb = !is.na(admits_bfb)), FUN = function(x) round(mean(x), 3)))
}
