# Harder simulation study: 22 chromosomes, measurement noise, smaller WGD
# clades, 10 replicates, and WGD-labelled reconstruction.
#   Rscript dev/wgd/run_experiments_hard.R [output.rds]
suppressMessages(devtools::load_all(".", quiet = TRUE))
Rcpp::sourceCpp("dev/wgd/wgd_proto.cpp")
source("dev/wgd/wgd_proto.R")
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "wgd_results_hard.rds"
chrs <- as.character(1:22)
n_reps <- 10
sim_args <- list(chromosomes = chrs, bfb_allele = "1:A", lambda = 1)

# ---- Measurement noise --------------------------------------------------------
# Per cell, chromosome and allele: each bin changes by +-1 with probability
# p_bin, and with probability p_seg a spurious segment of 1-5 bins changes
# by +-1. Copy numbers stay >= 0; CN = A + B.
add_noise <- function(data, p_bin, p_seg) {
  if (p_bin == 0 && p_seg == 0) return(data)
  dt <- as.data.table(data)
  setorder(dt, cell_id, chr, start)
  perturb <- function(v) {
    n <- length(v)
    flip <- runif(n) < p_bin
    v[flip] <- v[flip] + sample(c(-1, 1), sum(flip), replace = TRUE)
    if (runif(1) < p_seg) {
      len <- sample(1:5, 1); s <- sample(seq_len(max(1, n - len + 1)), 1)
      idx <- s:min(n, s + len - 1)
      v[idx] <- v[idx] + sample(c(-1, 1), 1)
    }
    pmax(v, 0)
  }
  dt[, `:=`(A = perturb(A), B = perturb(B)), by = .(cell_id, chr)]
  dt[, CN := A + B]
  dplyr::as_tibble(dt)
}

noise_levels <- list(none = c(0, 0), low = c(0.01, 0.1), high = c(0.03, 0.3))

# ---- Scenarios -----------------------------------------------------------------
from_sim <- function(sim) {
  tw <- true_wgd(sim)
  list(data = sim$cna_data, tree = sim$tree, n_wgd = tw$n_wgd, events = tw$events)
}
draw <- function(make, accept, start_seed) {
  out <- list(); seed <- start_seed
  while (length(out) < n_reps) {
    set.seed(seed)
    s <- tryCatch(suppressWarnings(make()), error = function(e) NULL)
    if (!is.null(s) && accept(s)) out[[length(out) + 1]] <- c(s, seed = seed)
    seed <- seed + 1
  }
  out
}
frac <- function(s, k) mean(s$n_wgd >= k)
base <- list(max_cells = 300, subsample = 150)

scenarios <- list(
  no_wgd = draw(function() from_sim(do.call(bridge_sim, c(sim_args, base, list(wgd_available = 0)))),
                function(s) TRUE, 1000),
  one_wgd = draw(function() from_sim(do.call(bridge_sim, c(sim_args, base, list(
                   wgd_available = 1, wgd_probability = 0.05)))),
                 function(s) frac(s, 1) >= 0.05 && frac(s, 1) <= 0.85, 2000),
  parallel_wgd = draw(function() from_sim(do.call(bridge_sim, c(sim_args, base, list(
                        wgd_available = 2, wgd_probability = 0.05)))),
                      function(s) length(s$events) == 2 && all(lengths(s$events) >= 4) &&
                        frac(s, 1) <= 0.85, 3000),
  nested_wgd = draw(function() {
      s1 <- do.call(bridge_sim, c(sim_args, list(max_cells = 8, wgd_available = 1,
                                                   wgd_probability = 0.5, return_phylo = FALSE)))
      k1 <- true_wgd(s1)$n_wgd
      s2 <- do.call(bridge_sim, c(sim_args, base, list(initial_sequences = s1$cells,
                                                        wgd_available = 1, wgd_probability = 0.05,
                                                        first_round_of_bfb = FALSE)))
      founder_k <- stats::setNames(k1, paste0("cell_", seq_along(k1)))
      tw <- true_wgd(s2)
      h <- s2$cell_history; par <- stats::setNames(h$parent_id, h$cell_id)
      founder_of <- vapply(names(tw$n_wgd), function(cell) {
        node <- cell; while (par[[node]] != "root") node <- par[[node]]; node }, "")
      list(data = s2$cna_data, tree = s2$tree, n_wgd = tw$n_wgd + founder_k[founder_of],
           events = tw$events)
    },
    function(s) all(c(0, 1, 2) %in% s$n_wgd) && sum(s$n_wgd == 2) >= 4, 4000)
)

# ---- Evaluation ------------------------------------------------------------------
evaluate <- function(s, data_noisy, alleles, mode) {
  data <- dplyr::bind_rows(data_noisy, create_diploid_data(data_noisy))
  t_dist <- system.time({
    X <- process_input_data(data, chrs, alleles, 0, 0)
    cls <- estimate_wgd_class(data)
    comp <- wgd_components(X, chrs, alleles)
  })[["elapsed"]]
  out <- list()
  for (m in mode) {
    d <- wgd_combine(comp, w = 1, mode = m, wgd_class = cls)
    tree <- tree_from_distance(d$D)
    cells <- names(s$n_wgd); k <- s$n_wgd[cells]
    cross <- outer(k, k, "!="); same <- outer(k, k, "==") & !diag(length(cells))
    D <- d$D[cells, cells]; W <- d$winner[cells, cells]
    groups <- c(s$events, if (any(k == 2)) list(k2 = cells[k == 2]))
    groups <- groups[lengths(groups) >= 2]

    # WGD-aware reconstruction and labels
    kn <- assign_wgd_counts(tree, cls[tree$tip.label])
    rec <- reconstruct_wgd(tree, X, chrs, alleles, kn)
    true_sets <- c(s$events, if (any(k == 2)) list(k2 = cells[k == 2]))
    sc <- score_wgd_labels(true_sets, wgd_clade_sets(tree, rec))
    # standard reconstruction on the same tree, for comparison
    std <- compute_reconstructions(list(tree = tree, all_input_Xs = X, g_dist_func = "greedy_fast",
                                        b_dist_func = "bfb_fast"), chrs, alleles)

    out[[m]] <- data.frame(
      mode = m,
      rf = phangorn::RF.dist(s$tree, tree, normalize = TRUE),
      wgd_clades_monophyletic = if (length(groups)) mean(sapply(groups, function(g) ape::is.monophyletic(tree, g))) else NA,
      mean_D_cross = if (any(cross)) mean(D[cross]) else NA,
      mean_D_within = mean(D[same]),
      wgd_wins_within = mean(W[same] > 0),
      class_accuracy = mean(cls[cells] == k),
      label_jaccard = if (length(true_sets)) mean(sc$best_jaccard) else NA,
      n_true_wgd = sc$n_true, n_pred_wgd = sc$n_pred,
      events_standard = sum(std$edge.length),
      events_wgd_aware = rec$total_events,
      distance_seconds = t_dist)
  }
  do.call(rbind, out)
}

results <- list()
for (sc in names(scenarios)) for (rep in seq_along(scenarios[[sc]])) {
  s <- scenarios[[sc]][[rep]]
  for (nl in names(noise_levels)) {
    set.seed(s$seed + 7)
    noisy <- add_noise(s$data, noise_levels[[nl]][1], noise_levels[[nl]][2])
    for (alleles in list(c("A", "B"), "CN")) {
      r <- evaluate(s, noisy, alleles, c("none", "all", "ploidy"))
      results[[length(results) + 1]] <- cbind(scenario = sc, rep = rep, seed = s$seed, noise = nl,
                                              data = paste(alleles, collapse = "+"), r)
    }
  }
  message(sc, " rep ", rep, " done")
  saveRDS(do.call(rbind, results), out_file)
}
results <- do.call(rbind, results)
saveRDS(results, out_file)
