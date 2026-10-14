# Do root-to-tip lengths of the fitted tree track the true number of events?
# In simulations the history records every event (bfb, wgd, amp, del) of every
# division, so each sampled cell's true event count is the sum over its
# lineage. Compared against:
#   depth  = root-to-tip path length of fit()$tree (reconstructed events per
#            edge, plus the WGD penalty when wgd = TRUE);
#   D_dipl = the distance from the cell to the synthetic diploid in fit()$D.
#   Rscript dev/branch_length/events_vs_depth.R [output.rds]
suppressMessages(library(bridges))
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "events_vs_depth.rds"
chrs <- as.character(1:10)
n_reps <- 5
types <- c("amp", "del", "bfb", "wgd")


scenarios <- list(
  cna_only = list(bfb_prob = 0,   first_round_of_bfb = FALSE, wgd_available = 0),
  bfb      = list(bfb_prob = 0.5, first_round_of_bfb = TRUE,  wgd_available = 0),
  wgd      = list(bfb_prob = 0,   first_round_of_bfb = FALSE, wgd_available = 1),
  bfb_wgd  = list(bfb_prob = 0.5, first_round_of_bfb = TRUE,  wgd_available = 1)
)

# True events per sampled cell, by type, summed along its lineage
true_events <- function(sim) {
  h <- sim$cell_history
  parent <- stats::setNames(h$parent_id, h$cell_id)
  ev <- strsplit(h$cn_event, ",", fixed = TRUE)
  per_row <- t(vapply(ev, function(e) vapply(types, function(t) sum(e == t), 0), numeric(length(types))))
  rownames(per_row) <- h$cell_id
  cells <- names(sim$cells)
  out <- matrix(0, length(cells), length(types), dimnames = list(cells, types))
  for (cell in cells) {
    node <- cell
    while (!is.na(node) && node != "root") {
      out[cell, ] <- out[cell, ] + per_row[node, ]
      node <- parent[[node]]
    }
  }
  out
}

draw <- function(sc, lambda, seed) {
  repeat {
    set.seed(seed)
    sim <- suppressWarnings(do.call(bridge_sim, c(list(chromosomes = chrs, bfb_allele = "1:A", lambda = lambda,
             max_cells = 300, subsample = 150, wgd_probability = 0.05), sc)))
    if (length(sim$cells) >= 100) {
      tr <- true_events(sim)
      # with WGD, keep runs where 15-85% of the cells are doubled
      f <- mean(tr[, "wgd"] > 0)
      if (sc$wgd_available == 0 || (f >= 0.15 && f <= 0.85)) return(list(sim = sim, truth = tr, seed = seed))
    }
    seed <- seed + 1
  }
}

rows <- list()
for (sc in names(scenarios)) for (lambda in c(0.5, 1, 2)) for (rep in seq_len(n_reps)) {
  s <- draw(scenarios[[sc]], lambda, 10000 * which(names(scenarios) == sc) + 1000 * lambda + 7 * rep)
  for (w in c(FALSE, TRUE)) {
    res <- fit(s$sim$cna_data, chromosomes = chrs, wgd = w, verbose = FALSE)
    cells <- res$tree$tip.label
    depth <- ape::node.depth.edgelength(res$tree)[seq_along(cells)]
    rows[[length(rows) + 1]] <- data.frame(
      scenario = sc, lambda = lambda, rep = rep, seed = s$seed, fit_wgd = w, cell = cells,
      s$truth[cells, , drop = FALSE], true_total = rowSums(s$truth[cells, , drop = FALSE]),
      depth = depth, D_dipl = res$D[cells, "diploid"], row.names = NULL)
  }
  message(sc, " lambda ", lambda, " rep ", rep, " done")
  saveRDS(do.call(rbind, rows), out_file)
}
