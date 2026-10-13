# Is BFB that happens AFTER a WGD still detected by fit(wgd = TRUE)?
# Two-stage simulations: stage 1 grows 8 founders without BFB, with a WGD
# (all founders doubled, or a mix); stage 2 starts from those founders and
# runs BFB on chromosome 1, allele A (or no BFB, as a control).
# Compares detect_bfb_branches() for wgd = FALSE, wgd = TRUE (current
# reconstruction, greedy target doubled inside doubled clades), and
# wgd = TRUE with the target not doubled.
#   Rscript dev/bfb_wgd/bfb_after_wgd_sim.R [output.rds]
suppressMessages(library(bridges))
source("dev/wgd/wgd_proto.R")   # true_wgd()
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "bfb_after_wgd_sim.rds"
chrs <- as.character(1:8)
n_reps <- 6

# Variant of reconstruct_tree() whose greedy-ancestor target is not doubled
lines <- readLines("R/bfb_mapping_and_detection.R")
s <- grep("^reconstruct_tree <- function", lines)
e <- s - 1 + which(lines[s:length(lines)] == "}")[1]
make <- function(txt) { f <- eval(parse(text = txt)); environment(f) <- asNamespace("bridges"); f }
rt_current <- make(lines[s:e])
rt_plain_target <- make(gsub("target_val * 2^k_here", "target_val", lines[s:e], fixed = TRUE))
use_reconstruction <- function(f) {
  ns <- asNamespace("bridges")
  unlockBinding("reconstruct_tree", ns); assign("reconstruct_tree", f, envir = ns); lockBinding("reconstruct_tree", ns)
}

scenarios <- list(
  all_doubled_bfb   = list(wgd_probability = 1,   bfb_prob = 0.5),
  mixed_bfb         = list(wgd_probability = 0.3, bfb_prob = 0.5),
  all_doubled_nobfb = list(wgd_probability = 1,   bfb_prob = 0),
  mixed_nobfb       = list(wgd_probability = 0.3, bfb_prob = 0)
)

simulate <- function(sc, seed) {
  repeat {
    set.seed(seed)
    s1 <- suppressWarnings(bridge_sim(chromosomes = chrs, bfb_allele = "1:A", max_cells = 8, lambda = 1,
                                      bfb_prob = 0, first_round_of_bfb = FALSE, wgd_available = 1,
                                      wgd_probability = sc$wgd_probability, return_phylo = FALSE))
    k1 <- true_wgd(s1)$n_wgd
    ok_founders <- length(s1$cells) == 8 && (if (sc$wgd_probability == 1) all(k1 == 1) else any(k1 == 1) && any(k1 == 0))
    if (ok_founders) {
      s2 <- suppressWarnings(bridge_sim(initial_sequences = s1$cells, chromosomes = chrs, bfb_allele = "1:A",
                                        hotspot = list(chr = "1:A", pos = 100), lambda = 1,
                                        bfb_prob = sc$bfb_prob, first_round_of_bfb = FALSE, wgd_available = 0,
                                        max_cells = 300, subsample = 150, return_phylo = FALSE))
      if (length(s2$cells) >= 100) return(list(data = s2$cna_data, seed = seed, frac_doubled = mean(k1 == 1)))
    }
    seed <- seed + 1
  }
}

rows <- list()
for (sc in names(scenarios)) for (rep in seq_len(n_reps)) {
  s <- simulate(scenarios[[sc]], 5000 * which(names(scenarios) == sc) + 37 * rep)
  for (alleles in list(c("A", "B"), "CN")) for (method in c("wgd = FALSE", "wgd = TRUE (current)", "wgd = TRUE, target not doubled")) {
    use_reconstruction(if (method == "wgd = TRUE, target not doubled") rt_plain_target else rt_current)
    res <- fit(s$data, chromosomes = chrs, alleles = alleles, wgd = method != "wgd = FALSE", verbose = FALSE)
    det <- detect_bfb_branches(res)
    truth <- det$chr == "1" & det$allele %in% c("A", "CN")
    rows[[length(rows) + 1]] <- data.frame(
      scenario = sc, rep = rep, data = paste(alleles, collapse = "+"), method = method,
      frac_founders_doubled = s$frac_doubled,
      bfb_1A_called = det$adj.pval[truth] < 0.05, bfb_1A_frac = det$mean[truth],
      false_calls = sum(det$adj.pval[!truth] < 0.05), false_frac_mean = mean(det$mean[!truth]))
  }
  use_reconstruction(rt_current)
  message(sc, " rep ", rep, " done")
  saveRDS(do.call(rbind, rows), out_file)
}
