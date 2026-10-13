# OV2295: ploidy-normalised distance (variant 1) and clone-level shared vs
# parallel WGD test (variant 3).
#   Rscript dev/wgd/real_alternatives.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
Rcpp::sourceCpp("dev/wgd/wgd_proto.cpp")
source("dev/wgd/wgd_proto.R")
source("dev/wgd/real_data_prep.R")
chrs <- as.character(1:22)

for (f in c(4, 1)) {
  rd <- load_ov2295(f)
  data <- dplyr::bind_rows(rd$data, create_diploid_data(rd$data))
  X <- process_input_data(data, chrs, "CN", 0, 0)
  cls <- estimate_wgd_class(data)
  clones <- rd$clones[names(rd$clones) %in% rownames(X[[1]]$CN)]
  cat(sprintf("\n==== bins x%d (%.1f Mb) ====\n", f, 0.5 * f))

  # Variant 1
  t_s0 <- system.time(S0 <- compute_fast_distances(X, chrs, "CN")$D)[["elapsed"]]
  t_norm <- system.time(norm <- wgd_normalised_distance(X, chrs, "CN", cls)$D)[["elapsed"]]
  tree <- tree_from_distance(norm)
  coh <- clone_coherence(tree, clones)
  kn <- assign_wgd_counts(tree, cls[tree$tip.label])
  rec <- reconstruct_wgd(tree, X, chrs, "CN", kn)
  sets <- wgd_clade_sets(tree, rec)
  big <- order(lengths(sets), decreasing = TRUE)
  nd <- names(clones)[clones %in% c("A", "B", "C", "D")]
  wg <- setdiff(names(clones), nd)
  cat(sprintf("normalised distance: %.1fs (S0 alone %.1fs) | mean clone best-F1 %.3f | mean D near-diploid vs others %.1f | WGD edges %d | events %.0f\n",
              t_norm, t_s0, weighted.mean(coh$best_f1, coh$n), mean(norm[nd, wg]), length(sets), rec$total_events))
  for (i in head(big, 3)) {
    tb <- table(clones[sets[[i]]])
    cat(sprintf("  WGD clade of %d cells: %s\n", length(sets[[i]]), paste(names(tb), tb, sep = ":", collapse = " ")))
  }

  # Variant 3: R = near-diploid cells; P, Q = doubled clones or groups
  R <- names(clones)[clones %in% c("A", "B", "C", "D") & cls[names(clones)] == 0]
  doubled <- function(cl) names(clones)[clones %in% cl & cls[names(clones)] == 1]
  tests <- list(c("E", "F"), c("G", "H"), c("G", "I"), c("H", "I"),
                c("E", "G"), c("E", "I"), c("F", "H"), list(c("E", "F"), c("G", "H", "I")))
  cat("clone-level test (delta = parallel - shared; > 0 favours one shared WGD):\n")
  for (t in tests) {
    P <- doubled(t[[1]]); Q <- doubled(t[[2]])
    r <- shared_vs_parallel(X, chrs, "CN", R, P, Q, w = 1)
    cat(sprintf("  %-4s vs %-4s  shared %.0f  parallel %.0f  delta %+.0f  -> %s\n",
                paste(t[[1]], collapse = ""), paste(t[[2]], collapse = ""),
                r[["shared"]], r[["parallel"]], r[["delta"]],
                ifelse(r[["delta"]] > 0, "one shared WGD", "two separate WGDs")))
  }
}
