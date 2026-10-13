# OV2295: event model (with BFB distance and exact WGD terms) vs the current
# greedy model with the full WGD search, at 2 Mb.
#   Rscript dev/wgd/real_event_model.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
Rcpp::sourceCpp("dev/wgd/wgd_proto.cpp")
Rcpp::sourceCpp("dev/wgd/event_model.cpp")
source("dev/wgd/wgd_proto.R")
source("dev/wgd/event_model.R")
source("dev/wgd/real_data_prep.R")
chrs <- as.character(1:22)

rd <- load_ov2295(4)
data <- dplyr::bind_rows(rd$data, create_diploid_data(rd$data))
X <- process_input_data(data, chrs, "CN", 0, 0)
cls <- estimate_wgd_class(data)
clones <- rd$clones[names(rd$clones) %in% rownames(X[[1]]$CN)]
nd <- names(clones)[clones %in% c("A", "B", "C", "D")]
wg <- setdiff(names(clones), nd)

t_g <- system.time(comp_g <- wgd_components(X, chrs, "CN"))[["elapsed"]]
t_e <- system.time(comp_e <- event_components(X, chrs, "CN", use_bfb = TRUE))[["elapsed"]]
cat(sprintf("component time at 2 Mb: greedy %.0fs, event + BFB %.0fs\n", t_g, t_e))

for (m in c("greedy_wgd", "event_bfb_no_wgd", "event_bfb_wgd")) {
  comp <- if (m == "greedy_wgd") comp_g else comp_e
  d <- wgd_combine(comp, 1, if (m == "event_bfb_no_wgd") "none" else "all", cls)
  tree <- tree_from_distance(d$D)
  coh <- clone_coherence(tree, clones)
  kn <- assign_wgd_counts(tree, cls[tree$tip.label])
  rec <- if (m == "greedy_wgd") reconstruct_wgd(tree, X, chrs, "CN", kn) else
    reconstruct_wgd_event(tree, X, chrs, "CN", kn)
  sets <- wgd_clade_sets(tree, rec)
  big <- order(lengths(sets), decreasing = TRUE)
  cat(sprintf("\n-- %s --\nmean clone best-F1 %.3f | D near-diploid vs others %.1f (within near-diploid %.1f) | WGD edges %d\n",
              m, weighted.mean(coh$best_f1, coh$n), mean(d$D[nd, wg]),
              mean(d$D[nd, nd][upper.tri(d$D[nd, nd])]), length(sets)))
  print(data.frame(clone = coh$clone, n = coh$n, best_f1 = round(coh$best_f1, 3)), row.names = FALSE)
  for (i in head(big, 3)) {
    tb <- table(clones[sets[[i]]])
    cat(sprintf("  WGD clade of %d cells: %s\n", length(sets[[i]]), paste(names(tb), tb, sep = ":", collapse = " ")))
  }
}
