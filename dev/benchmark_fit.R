# Time each stage of fit() on simulated data. Run from the package root:
#   Rscript dev/benchmark_fit.R
# Uses internal functions, so it loads the package with devtools::load_all().
suppressMessages(devtools::load_all(".", quiet = TRUE))

tm <- function(expr) round(system.time(expr)[["elapsed"]], 2)

configs <- list(
  list(cells = 400, chrs = as.character(1:6)),
  list(cells = 1000, chrs = as.character(1:22))
)

for (cfg in configs) {
  set.seed(7)
  sim <- suppressWarnings(bridge_sim(
    chromosomes = cfg$chrs, bfb_allele = "1:A", max_cells = 2 * cfg$cells,
    subsample = cfg$cells, lambda = 1, return_phylo = FALSE
  ))
  d <- sim$cna_data
  data <- dplyr::bind_rows(d, create_diploid_data(d))

  t_pre <- tm(X <- process_input_data(data, cfg$chrs, c("A", "B"), 0, 0))
  t_dist <- tm(dists <- compute_fast_distances(X, cfg$chrs, c("A", "B")))
  D <- dists$D
  diag(D) <- 0
  t_nj <- tm(tree <- ape::nj(D))
  tree <- ape::drop.tip(ape::root(tree, "diploid", resolve.root = TRUE), "diploid")
  partial_fit <- list(tree = tree, all_input_Xs = X,
                      g_dist_func = "greedy_fast", b_dist_func = "bfb_fast")
  t_rec <- tm(compute_reconstructions(partial_fit, cfg$chrs, c("A", "B")))
  t_fit <- tm(res <- fit(d, chromosomes = cfg$chrs, verbose = FALSE))

  cat(sprintf(
    "%d cells x %d chromosomes: fit() %.2f s = pre-processing %.2f | distances %.2f | nj %.2f | reconstruction %.2f | result %.0f MB\n",
    cfg$cells, length(cfg$chrs), t_fit, t_pre, t_dist, t_nj, t_rec,
    as.numeric(object.size(res)) / 2^20
  ))
}
