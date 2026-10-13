# Calibration and power of detect_bfb_clades() vs detect_bfb_branches().
# BFB (when present) is only on chromosome 1, allele A.
#   Rscript dev/bfb_clades/clade_test_sim.R [output.rds]
# Runs on the installed package.
suppressMessages(library(bridges))
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "clade_test_sim.rds"
chrs <- as.character(1:8)
n_reps <- 10

sim1 <- function(...) suppressWarnings(bridge_sim(chromosomes = chrs, bfb_allele = "1:A",
                                                    hotspot = list(chr = "1:A", pos = 100), lambda = 1, ...))

# Returns list(data, bfb_cells) or NULL (e.g. extinction)
make <- function(scenario, seed) {
  set.seed(seed)
  if (scenario %in% c("no_bfb", "no_bfb_with_wgd", "whole_population_bfb")) {
    s <- sim1(max_cells = 300, subsample = 150,
              bfb_prob = if (scenario == "whole_population_bfb") 0.5 else 0,
              first_round_of_bfb = scenario == "whole_population_bfb",
              wgd_available = if (scenario == "no_bfb_with_wgd") 1 else 0, wgd_probability = 0.05,
              return_phylo = FALSE)
    if (length(s$cells) < 100) return(NULL)
    return(list(data = s$cna_data, bfb_cells = if (scenario == "whole_population_bfb") unique(s$cna_data$cell_id) else character(0)))
  }
  # subpopulation BFB: 8 founders without BFB; 2 of them then evolve with BFB,
  # the other 6 without; the two populations are sequenced together
  s1 <- sim1(max_cells = 8, bfb_prob = 0, first_round_of_bfb = FALSE, return_phylo = FALSE)
  if (length(s1$cells) != 8) return(NULL)
  a <- sim1(initial_sequences = s1$cells[1:2], max_cells = 60, subsample = 30, bfb_prob = 0.5,
            first_round_of_bfb = FALSE, return_phylo = FALSE)
  b <- sim1(initial_sequences = s1$cells[3:8], max_cells = 300, subsample = 120, bfb_prob = 0,
            first_round_of_bfb = FALSE, return_phylo = FALSE)
  if (length(a$cells) < 20 || length(b$cells) < 80) return(NULL)
  a$cna_data$cell_id <- paste0("bfb_", a$cna_data$cell_id)
  b$cna_data$cell_id <- paste0("other_", b$cna_data$cell_id)
  list(data = dplyr::bind_rows(a$cna_data, b$cna_data), bfb_cells = unique(a$cna_data$cell_id))
}

scenarios <- c("no_bfb", "no_bfb_with_wgd", "whole_population_bfb", "subpopulation_bfb")
rows <- list()
for (sc in scenarios) for (rep in seq_len(n_reps)) {
  seed <- 7000 * which(scenarios == sc) + 13 * rep
  repeat { s <- make(sc, seed); if (!is.null(s)) break; seed <- seed + 1 }
  res <- fit(s$data, chromosomes = chrs, wgd = sc == "no_bfb_with_wgd", verbose = FALSE)
  g <- detect_bfb_branches(res)
  t_cl <- system.time(cl <- detect_bfb_clades(res, n_sim = 999, seed = rep))[["elapsed"]]
  key_g <- paste(g$chr, g$allele); key_c <- paste(cl$chr, cl$allele)
  truth <- "1 A"
  j <- if (length(s$bfb_cells)) {
    best <- cl$clade_cells[[which(key_c == truth)]]
    length(intersect(best, s$bfb_cells)) / length(union(best, s$bfb_cells))
  } else NA
  rows[[length(rows) + 1]] <- data.frame(
    scenario = sc, rep = rep, seed = seed, n_cells = length(unique(s$data$cell_id)),
    n_bfb_cells = length(s$bfb_cells),
    global_1A = g$adj.pval[key_g == truth] < 0.05,
    clade_1A = cl$adj.pval[key_c == truth] < 0.05,
    global_false = sum(g$adj.pval[key_g != truth] < 0.05),
    clade_false = sum(cl$adj.pval[key_c != truth] < 0.05),
    n_tests_other = sum(key_c != truth),
    clade_jaccard_1A = j,
    clade_seconds = t_cl)
  message(sc, " rep ", rep, " done")
  saveRDS(do.call(rbind, rows), out_file)
}
