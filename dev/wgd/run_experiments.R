# Simulation study for the prototype WGD-aware distance.
#   Rscript dev/wgd/run_experiments.R [output.rds]
suppressMessages(devtools::load_all(".", quiet = TRUE))
Rcpp::sourceCpp("dev/wgd/wgd_proto.cpp")
source("dev/wgd/wgd_proto.R")

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "wgd_results.rds"
chrs <- as.character(1:10)
n_reps <- 6

sim_args <- list(chromosomes = chrs, bfb_allele = "1:A", lambda = 1)

# ---- Scenario generators: return list(data, tree, n_wgd, events) ----------

from_sim <- function(sim) {
  tw <- true_wgd(sim)
  list(data = sim$cna_data, tree = sim$tree, n_wgd = tw$n_wgd, events = tw$events)
}

# Keep drawing seeds until `accept()` holds; returns n_reps scenarios.
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

scenarios <- list(
  no_wgd = draw(function() from_sim(do.call(bridge_sim, c(sim_args, list(
      max_cells = 600, subsample = 150, wgd_available = 0)))),
    function(s) TRUE, 100),
  one_wgd = draw(function() from_sim(do.call(bridge_sim, c(sim_args, list(
      max_cells = 600, subsample = 150, wgd_available = 1, wgd_probability = 0.05)))),
    function(s) frac(s, 1) >= 0.15 && frac(s, 1) <= 0.85, 200),
  parallel_wgd = draw(function() from_sim(do.call(bridge_sim, c(sim_args, list(
      max_cells = 600, subsample = 150, wgd_available = 2, wgd_probability = 0.05)))),
    function(s) length(s$events) == 2 && all(lengths(s$events) >= 8) && frac(s, 1) <= 0.85, 300),
  nested_wgd = draw(function() {
      # Stage 1: 8 founders, some doubled once. Stage 2: one more WGD.
      s1 <- do.call(bridge_sim, c(sim_args, list(max_cells = 8, wgd_available = 1,
                                                   wgd_probability = 0.5, return_phylo = FALSE)))
      k1 <- true_wgd(s1)$n_wgd
      s2 <- do.call(bridge_sim, c(sim_args, list(initial_sequences = s1$cells, max_cells = 600,
                                                   subsample = 150, wgd_available = 1,
                                                   wgd_probability = 0.05, first_round_of_bfb = FALSE)))
      # initialize_from_sequences() renames founders cell_1..cell_n in order
      founder_k <- stats::setNames(k1, paste0("cell_", seq_along(k1)))
      tw <- true_wgd(s2)
      founder_of <- sapply(names(tw$n_wgd), function(cell) {
        h <- s2$cell_history; node <- cell
        while (h$parent_id[h$cell_id == node] != "root") node <- h$parent_id[h$cell_id == node]
        node
      })
      list(data = s2$cna_data, tree = s2$tree,
           n_wgd = tw$n_wgd + founder_k[founder_of], events = tw$events)
    },
    function(s) all(c(0, 1, 2) %in% s$n_wgd) && sum(s$n_wgd == 2) >= 8, 400)
)

# ---- Evaluation -----------------------------------------------------------

prepare <- function(s, alleles) {
  data <- dplyr::bind_rows(s$data, create_diploid_data(s$data))
  X <- process_input_data(data, chrs, alleles, 0, 0)
  list(comp = wgd_components(X, chrs, alleles), cls = estimate_wgd_class(data))
}

evaluate <- function(s, prep, mode, w) {
  cls <- prep$cls
  d <- wgd_combine(prep$comp, w = w, mode = mode, wgd_class = cls)
  tree <- tree_from_distance(d$D)

  cells <- names(s$n_wgd)
  k <- s$n_wgd[cells]
  cross <- outer(k, k, "!=")
  same <- outer(k, k, "==") & !diag(length(cells))
  D <- d$D[cells, cells]; W <- d$winner[cells, cells]
  # monophyly of every WGD clade with at least 2 sampled cells
  groups <- c(s$events, if (any(k == 2)) list(k2 = cells[k == 2]))
  groups <- groups[lengths(groups) >= 2]
  mono <- if (length(groups)) mean(sapply(groups, function(g) ape::is.monophyletic(tree, g))) else NA

  data.frame(
    rf = phangorn::RF.dist(s$tree, tree, normalize = TRUE),
    wgd_clades_monophyletic = mono,
    mean_D_cross = if (any(cross)) mean(D[cross]) else NA,
    mean_D_within = mean(D[same]),
    wgd_wins_cross = if (any(cross)) mean(W[cross] > 0) else NA,
    wgd_wins_within = mean(W[same] > 0),
    class_accuracy = mean(cls[cells] == k)
  )
}

configs <- rbind(
  data.frame(mode = "none",   w = NA),
  data.frame(mode = "all",    w = c(1, 3, 10)),
  data.frame(mode = "ploidy", w = c(1, 3, 10))
)

results <- list()
for (sc in names(scenarios)) for (rep in seq_along(scenarios[[sc]])) {
  s <- scenarios[[sc]][[rep]]
  for (alleles in list(c("A", "B"), "CN")) {
    t_prep <- system.time(prep <- prepare(s, alleles))[["elapsed"]]
    for (i in seq_len(nrow(configs))) {
      cfg <- configs[i, ]
      r <- evaluate(s, prep, cfg$mode, ifelse(is.na(cfg$w), 1, cfg$w))
      results[[length(results) + 1]] <- cbind(
        scenario = sc, rep = rep, seed = s$seed, data = paste(alleles, collapse = "+"),
        mode = cfg$mode, w = cfg$w, r, distance_seconds = t_prep)
    }
  }
  message(sc, " rep ", rep, " done")
}
results <- do.call(rbind, results)
saveRDS(results, out_file)

summary_tbl <- aggregate(
  cbind(rf, wgd_clades_monophyletic, mean_D_cross, mean_D_within, wgd_wins_cross,
        wgd_wins_within, distance_seconds) ~ scenario + data + mode + w,
  data = transform(results, w = ifelse(is.na(w), "-", w)),
  FUN = function(x) round(mean(x), 3), na.action = na.pass)
summary_tbl <- summary_tbl[order(summary_tbl$scenario, summary_tbl$data, summary_tbl$mode, summary_tbl$w), ]
print(summary_tbl, row.names = FALSE)
