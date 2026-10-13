# Benchmark of the alternative WGD-aware approaches against the full search:
#   variant 1: ploidy-normalised distance (wgd_normalised_distance())
#   variant 3: clone-level shared vs parallel WGD test (shared_vs_parallel())
# Same scenarios and seeds as run_experiments_hard.R, noise "none" and "low".
#   Rscript dev/wgd/run_alternatives.R [output.rds]
src <- readLines("dev/wgd/run_experiments_hard.R")
# reuse setup, noise model and scenario generation (everything before "Evaluation")
eval(parse(text = src[1:(grep("# ---- Evaluation", src, fixed = TRUE) - 1)]))
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "wgd_alternatives.rds"

score_tree <- function(s, D, X, alleles, cls) {
  tree <- tree_from_distance(D)
  cells <- names(s$n_wgd); k <- s$n_wgd[cells]
  kn <- assign_wgd_counts(tree, cls[tree$tip.label])
  rec <- reconstruct_wgd(tree, X, chrs, alleles, kn)
  true_sets <- c(s$events, if (any(k == 2)) list(k2 = cells[k == 2]))
  sc <- score_wgd_labels(true_sets, wgd_clade_sets(tree, rec))
  cross <- outer(k, k, "!=")
  data.frame(rf = phangorn::RF.dist(s$tree, tree, normalize = TRUE),
             label_jaccard = if (length(true_sets)) mean(sc$best_jaccard) else NA,
             n_true_wgd = sc$n_true, n_pred_wgd = sc$n_pred,
             mean_D_cross = if (any(cross)) mean(D[cells, cells][cross]) else NA,
             events_wgd_aware = rec$total_events)
}

# Groups for the clone-level test, with the true answer
test_groups <- function(s, scenario) {
  cells <- names(s$n_wgd); k <- s$n_wgd[cells]
  R <- cells[k == 0]
  if (scenario == "one_wgd") {
    wgd_cells <- cells[k == 1]
    tree <- s$tree
    mrca <- ape::getMRCA(tree, wgd_cells)
    kids <- tree$edge[tree$edge[, 1] == mrca, 2]
    tips <- function(node) if (node <= length(tree$tip.label)) tree$tip.label[node] else
      ape::extract.clade(tree, node)$tip.label
    sets <- lapply(kids, function(x) intersect(tips(x), wgd_cells))
    sets <- sets[order(lengths(sets), decreasing = TRUE)][1:2]
    list(R = R, P = sets[[1]], Q = sets[[2]], truth = "shared")
  } else if (scenario == "parallel_wgd") {
    ev <- s$events[order(lengths(s$events), decreasing = TRUE)][1:2]
    list(R = R, P = ev[[1]], Q = ev[[2]], truth = "parallel")
  } else NULL
}

results <- list(); tests <- list()
for (sc in c("no_wgd", "one_wgd", "parallel_wgd", "nested_wgd")) for (rep in seq_along(scenarios[[sc]])) {
  s <- scenarios[[sc]][[rep]]
  for (nl in c("none", "low")) {
    set.seed(s$seed + 7)
    noisy <- add_noise(s$data, noise_levels[[nl]][1], noise_levels[[nl]][2])
    for (alleles in list(c("A", "B"), "CN")) {
      data <- dplyr::bind_rows(noisy, create_diploid_data(noisy))
      X <- process_input_data(data, chrs, alleles, 0, 0)
      cls <- estimate_wgd_class(data)
      t_s0 <- system.time(S0 <- compute_fast_distances(X, chrs, alleles)$D)[["elapsed"]]
      diag(S0) <- 0
      t_full <- system.time(full <- wgd_combine(wgd_components(X, chrs, alleles), 1, "all", cls)$D)[["elapsed"]]
      t_norm <- system.time(norm <- wgd_normalised_distance(X, chrs, alleles, cls, w = 1)$D)[["elapsed"]]
      for (m in c("none", "full", "normalised")) {
        D <- switch(m, none = S0, full = full, normalised = norm)
        results[[length(results) + 1]] <- cbind(
          scenario = sc, rep = rep, noise = nl, data = paste(alleles, collapse = "+"), method = m,
          score_tree(s, D, X, alleles, cls),
          distance_seconds = switch(m, none = t_s0, full = t_full, normalised = t_norm))
      }
      g <- test_groups(s, sc)
      if (!is.null(g) && length(g$P) && length(g$Q) && length(g$R)) {
        r <- shared_vs_parallel(X, chrs, alleles, g$R, g$P, g$Q, w = 1)
        tests[[length(tests) + 1]] <- data.frame(
          scenario = sc, rep = rep, noise = nl, data = paste(alleles, collapse = "+"),
          size_P = length(g$P), size_Q = length(g$Q), truth = g$truth,
          shared = r[["shared"]], parallel = r[["parallel"]], delta = r[["delta"]],
          call = ifelse(r[["delta"]] > 0, "shared", ifelse(r[["delta"]] < 0, "parallel", "tie")))
      }
    }
  }
  message(sc, " rep ", rep, " done")
  saveRDS(list(trees = do.call(rbind, results), tests = do.call(rbind, tests)), out_file)
}
