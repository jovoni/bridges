# How does WGD affect BFB detection? Simulations with known truth: BFB only
# on chromosome 1, allele A.
#   Rscript dev/bfb_wgd/bfb_vs_wgd_sim.R [output.rds]
# Runs on the installed package.
suppressMessages(library(bridges))
source("dev/wgd/wgd_proto.R")   # true_wgd()

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "bfb_vs_wgd_sim.rds"
chrs <- as.character(1:10)
n_reps <- 5

scenarios <- list(
  bfb_and_wgd = list(bfb_prob = 0.5, first_round_of_bfb = TRUE, wgd_available = 1),
  wgd_only    = list(bfb_prob = 0,   first_round_of_bfb = FALSE, wgd_available = 1),
  bfb_only    = list(bfb_prob = 0.5, first_round_of_bfb = TRUE, wgd_available = 0)
)

# First simulation from start_seed onwards that did not go extinct and (with
# WGD) has 15-85% doubled cells
draw <- function(sc, start_seed) {
  seed <- start_seed
  repeat {
    set.seed(seed)
    sim <- suppressWarnings(do.call(bridge_sim, c(list(chromosomes = chrs, bfb_allele = "1:A", lambda = 1,
             max_cells = 300, subsample = 150, wgd_probability = 0.05), sc)))
    if (length(sim$cells) > 0) {
      f <- mean(true_wgd(sim)$n_wgd > 0)
      if (sc$wgd_available == 0 || (f >= 0.15 && f <= 0.85)) return(list(sim = sim, seed = seed))
    }
    seed <- seed + 1
  }
}

tree_rows <- list(); cell_rows <- list()
for (sc in names(scenarios)) {
  for (rep in seq_len(n_reps)) {
    s <- draw(scenarios[[sc]], 1000 * which(names(scenarios) == sc) + 10 * rep)
    sim <- s$sim
    doubled <- names(which(true_wgd(sim)$n_wgd > 0))

    # Tree-based detection, with and without the WGD option
    for (w in c(FALSE, TRUE)) {
      res <- fit(sim$cna_data, chromosomes = chrs, wgd = w, verbose = FALSE)
      det <- detect_bfb_branches(res)
      det$truth <- det$chr == "1" & det$allele == "A"
      tree_rows[[length(tree_rows) + 1]] <- data.frame(
        scenario = sc, rep = rep, wgd = w,
        true_bfb_called = if (any(det$truth)) det$adj.pval[det$truth] < 0.05 else NA,
        true_bfb_frac = if (any(det$truth)) det$mean[det$truth] else NA,
        false_calls = sum(det$adj.pval[!det$truth] < 0.05),
        false_frac_mean = mean(det$mean[!det$truth]),
        n_tests = nrow(det))
    }

    # Per-cell detection: raw, and with doubled cells halved first
    halve <- sim$cna_data
    est <- bridges:::estimate_wgd_counts(sim$cna_data)
    hi <- halve$cell_id %in% names(est)[est > 0]
    halve$A[hi] <- floor(halve$A[hi] / 2); halve$B[hi] <- floor(halve$B[hi] / 2)
    halve$CN <- halve$A + halve$B
    for (variant in c("raw", "halved")) {
      d <- if (variant == "raw") sim$cna_data else halve
      det <- bfb_detect_batch(d, alleles = c("A", "B"), genome = "hg19")
      det <- det[!is.na(det$admits_bfb), ]
      det$bfb_truth <- det$chr == "1" & det$allele == "A"
      det$doubled <- det$cell_id %in% doubled
      agg <- aggregate(admits_bfb ~ bfb_truth + doubled, data = det, FUN = mean)
      n <- aggregate(admits_bfb ~ bfb_truth + doubled, data = det, FUN = length)
      cell_rows[[length(cell_rows) + 1]] <- cbind(scenario = sc, rep = rep, variant = variant,
                                                  agg, n_informative = n$admits_bfb)
    }
    message(sc, " rep ", rep, " done")
    saveRDS(list(tree = do.call(rbind, tree_rows), cells = do.call(rbind, cell_rows)), out_file)
  }
}
