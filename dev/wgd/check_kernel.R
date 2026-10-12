# Checks for the prototype WGD kernel: known cases and exhaustive search.
#   Rscript dev/wgd/check_kernel.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
Rcpp::sourceCpp("dev/wgd/wgd_proto.cpp")
G <- function(a, b) greedy_distance_cpp(as.numeric(a), as.numeric(b), 1)$cost
x <- c(1, 1, 3, 3, 1, 0, 0, 2, 2, 1)
chk <- function(label, got, want) cat(sprintf("%-45s got %s, expected %s  %s\n", label, got, want, if (got == want) "OK" else "MISMATCH"))
chk("W(x -> 2x)", wgd_cost(x, 2 * x, 2), 0)
chk("W(x -> 2x, one focal change)", wgd_cost(x, replace(2 * x, 3, 5), 2), 1)
# 2x followed by halving the tail 4-4-2 -> 2-2-1: not one greedy event; exact minimum is 2
chk("W(x -> 2x, tail halved)", wgd_cost(x, c(2 * x[1:5], x[6:10]), 2), 2)
chk("W(x -> 4x), f = 4", wgd_cost(x, 4 * x, 4), 0)
chk("W(flat 1 -> flat 2)", wgd_cost(rep(1, 5), rep(2, 5), 2), 0)
y <- replace(x, 3:4, 1)   # x with a focal loss before WGD on the other lineage
chk("W(y -> 2x), x one focal gain from y", wgd_cost(y, 2 * x, 2), 1)
cat("for reference, G(x, 2x) =", G(x, 2 * x), "\n\n")

# Exhaustive search over x on short vectors
set.seed(1)
gap <- c(); n_cases <- 400
for (k in seq_len(n_cases)) {
  n <- sample(2:4, 1); f <- sample(c(2, 4), 1)
  a <- sample(0:3, n, TRUE); b <- sample(0:9, n, TRUE)
  grid <- as.matrix(expand.grid(rep(list(0:5), n)))
  exact <- min(apply(grid, 1, function(xx) G(a, xx) + G(f * xx, b)))
  approx <- wgd_cost(a, b, f)
  gap <- c(gap, approx - exact)
}
cat(sprintf("exhaustive check, %d random short cases: approx == exact in %.1f%%, never below exact: %s, max gap %d, mean gap %.3f\n",
            n_cases, 100 * mean(gap == 0), all(gap >= 0), max(gap), mean(gap)))

# With no WGD term, the prototype distance must equal fit()'s distance matrix
source("dev/wgd/wgd_proto.R")
set.seed(3)
sim <- suppressWarnings(bridge_sim(chromosomes = c("1", "2"), bfb_allele = "1:A", max_cells = 120,
                                   subsample = 60, lambda = 1, wgd_available = 1,
                                   wgd_probability = 0.1, return_phylo = FALSE))
res <- fit(sim$cna_data, chromosomes = c("1", "2"), verbose = FALSE)
comp <- wgd_components(res$all_input_Xs, c("1", "2"), c("A", "B"))
D0 <- wgd_combine(comp, mode = "none")$D
cat("prototype mode = 'none' equals fit()$D:", identical(D0[rownames(res$D), colnames(res$D)], res$D), "\n")
