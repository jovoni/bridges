# Speed of fit(wgd = TRUE): WGD costs for all pairs vs only between cells
# whose WGD counts differ (default), and whether results are identical.
#   Rscript dev/benchmark_wgd.R
# Uses the OV2295 sample in dev/data if present.
# Runs on the installed package: devtools::load_all() compiles the C++ code
# without optimisation (-O0), which makes it several times slower. Install
# first, e.g. devtools::install().
suppressMessages(library(bridges))
library(data.table)

# CPU time (user + system) of this process: unaffected by other jobs or sleep
tm <- function(expr) { t <- system.time(expr); round(t[["user.self"]] + t[["sys.self"]], 2) }

# Measurement noise as in dev/wgd/run_experiments_hard.R
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

compare <- function(label, data, chrs, alleles) {
  t_none <- tm(r0 <- fit(data, chromosomes = chrs, alleles = alleles, verbose = FALSE))
  t_all <- tm(rA <- fit(data, chromosomes = chrs, alleles = alleles, wgd = TRUE,
                        wgd_all_pairs = TRUE, verbose = FALSE))
  t_fast <- tm(rF <- fit(data, chromosomes = chrs, alleles = alleles, wgd = TRUE, verbose = FALSE))
  same <- identical(rF$D, rA$D) && isTRUE(all.equal(rF$tree, rA$tree)) &&
    identical(rF$wgd$edges$cells, rA$wgd$edges$cells)
  cat(sprintf("%-38s fit() CPU: no WGD %6.1fs | WGD all pairs %6.1fs | WGD restricted %6.1fs | speed-up %4.1fx | identical: %s | WGD edges %d\n",
              label, t_none, t_all, t_fast, t_all / t_fast, same, nrow(rF$wgd$edges)))
}

chrs <- as.character(1:22)
for (cfg in list(list("sim, 1 WGD, no noise", 1, 0, 0, 201), list("sim, 1 WGD, low noise", 1, 0.01, 0.1, 201),
                 list("sim, 2 WGDs, no noise", 2, 0, 0, 305), list("sim, 2 WGDs, low noise", 2, 0.01, 0.1, 305))) {
  set.seed(cfg[[5]])
  sim <- suppressWarnings(bridge_sim(chromosomes = chrs, bfb_allele = "1:A", lambda = 1, max_cells = 300,
                                     subsample = 150, wgd_available = cfg[[2]], wgd_probability = 0.05,
                                     return_phylo = FALSE))
  d <- add_noise(sim$cna_data, cfg[[3]], cfg[[4]])
  compare(paste0(cfg[[1]], ", A+B"), d, chrs, c("A", "B"))
}

if (file.exists("dev/data/ov2295_cell_cn.csv.gz")) {
  source("dev/wgd/real_data_prep.R")
  for (f in c(4, 1)) {
    rd <- load_ov2295(f)
    compare(sprintf("OV2295, %.1f Mb, CN", 0.5 * f), rd$data, chrs, "CN")
  }
}
