G <- function(a, b) bridges:::greedy_distance_cpp(as.numeric(a), as.numeric(b), 1)$cost
W <- function(a, b, f) bridges:::wgd_cost_cpp(as.numeric(a), as.numeric(b), f)

test_that("wgd_cost_cpp gives the expected cost on known cases", {
  x <- c(1, 1, 3, 3, 1, 0, 0, 2, 2, 1)
  expect_equal(W(x, 2 * x, 2), 0)                    # pure doubling
  expect_equal(W(x, replace(2 * x, 3, 5), 2), 1)     # doubling + one focal change
  expect_equal(W(x, 4 * x, 4), 0)                    # two doublings
  expect_equal(W(rep(1, 5), rep(2, 5), 2), 0)
  expect_gt(G(x, 2 * x), 0)                          # without WGD this costs several events
})

test_that("wgd_cost_cpp is never below the exact minimum and usually equals it", {
  set.seed(31)
  gap <- replicate(300, {
    n <- sample(2:4, 1); f <- sample(c(2, 4), 1)
    a <- sample(0:3, n, TRUE); b <- sample(0:9, n, TRUE)
    grid <- as.matrix(expand.grid(rep(list(0:5), n)))
    exact <- min(apply(grid, 1, function(x) G(a, x) + G(f * x, b)))
    W(a, b, f) - exact
  })
  expect_true(all(gap >= 0))
  expect_gt(mean(gap == 0), 0.9)
  expect_lte(max(gap), 1)
})

test_that("wgd_ancestor_cpp returns an ancestor consistent with its cost", {
  set.seed(32)
  for (k in 1:200) {
    n <- sample(1:30, 1); f <- sample(c(2, 4), 1)
    a <- sample(0:4, n, TRUE); b <- sample(0:12, n, TRUE)
    r <- bridges:::wgd_ancestor_cpp(a, b, f)
    expect_equal(r$cost, W(a, b, f))
    expect_equal(r$after, G(f * r$x, b))
    expect_lte(G(a, r$x), r$steps)
  }
})

test_that("wgd_matrices_cpp matches per-pair wgd_cost_cpp, with repeated profiles", {
  set.seed(33)
  profiles <- matrix(sample(0:6, 6 * 15, TRUE), 6)
  X <- profiles[sample(6, 20, TRUE), ]
  counts <- sample(0:2, 20, TRUE)
  all <- bridges:::wgd_matrices_cpp(X, counts, all_pairs = TRUE)
  some <- bridges:::wgd_matrices_cpp(X, counts, all_pairs = FALSE)
  for (i in 1:20) for (j in 1:20) {
    if (i == j) next
    expect_equal(all$W1[i, j], W(X[i, ], X[j, ], 2))
    expect_equal(all$W2[i, j], W(X[i, ], X[j, ], 4))
    dk <- counts[j] - counts[i]
    expect_equal(some$W1[i, j], if (dk >= 1) all$W1[i, j] else Inf)
    expect_equal(some$W2[i, j], if (dk == 2) all$W2[i, j] else Inf)
  }
  expect_true(all(is.infinite(diag(all$W1))))
  expect_lt(some$n_computed, all$n_computed)
})

test_that("wgd_matrices_cpp keeps identical profiles with different WGD counts apart", {
  X <- rbind(c(1, 2, 1), c(1, 2, 1), c(2, 4, 2))
  m <- bridges:::wgd_matrices_cpp(X, c(0L, 1L, 1L))
  expect_equal(m$n_unique, 3)
  expect_equal(m$W1[1, 3], 0)               # third row is the first doubled
  expect_equal(m$W1[1, 2], W(X[1, ], X[2, ], 2))
})

test_that("assign_wgd_counts places WGDs on the expected edges", {
  shared <- ape::read.tree(text = "((d1,d2),(w1,w2));")
  k <- bridges:::assign_wgd_counts(shared, c(d1 = 0, d2 = 0, w1 = 1, w2 = 1))
  e <- bridges:::wgd_edges(shared, k)
  expect_equal(nrow(e), 1)
  expect_setequal(e$cells[[1]], c("w1", "w2"))

  parallel <- ape::read.tree(text = "((w1,d1),(w2,d2));")
  e <- bridges:::wgd_edges(parallel, bridges:::assign_wgd_counts(parallel, c(d1 = 0, d2 = 0, w1 = 1, w2 = 1)))
  expect_equal(nrow(e), 2)
  expect_setequal(unlist(e$cells), c("w1", "w2"))

  all_doubled <- bridges:::assign_wgd_counts(shared, c(d1 = 1, d2 = 1, w1 = 1, w2 = 1))
  expect_equal(all_doubled[5], 1L)   # root: WGD above the root
  expect_equal(nrow(bridges:::wgd_edges(shared, all_doubled)), 0)
})

wgd_sim <- function(seed, wgd_available) {
  set.seed(seed)
  suppressWarnings(bridge_sim(chromosomes = as.character(1:4), bfb_allele = "1:A", lambda = 1,
                              max_cells = 120, subsample = 50, wgd_available = wgd_available,
                              wgd_probability = 0.1))
}
doubled_cells <- function(sim) {
  h <- sim$cell_history
  names(sim$cells)[vapply(names(sim$cells), function(cell) {
    node <- cell; hit <- FALSE
    while (node != "root") {
      hit <- hit || h$wgd_event[h$cell_id == node]
      node <- h$parent_id[h$cell_id == node]
    }
    hit
  }, logical(1))]
}

test_that("fit(wgd = TRUE) finds a single WGD and its clade", {
  seed <- 1
  repeat {   # a simulation where the WGD clade is neither tiny nor everything
    sim <- wgd_sim(seed, 1); truth <- doubled_cells(sim)
    if (length(truth) >= 5 && length(truth) <= 40) break
    seed <- seed + 1
  }
  chrs <- as.character(1:4)
  r0 <- fit(sim$cna_data, chromosomes = chrs, verbose = FALSE)
  r1 <- fit(sim$cna_data, chromosomes = chrs, wgd = TRUE, verbose = FALSE)

  expect_null(r0$wgd)
  expect_equal(nrow(r1$wgd$edges), 1)
  expect_setequal(r1$wgd$edges$cells[[1]], truth)
  others <- setdiff(names(sim$cells), truth)
  expect_lt(mean(r1$D[truth, others]), mean(r0$D[truth, others]))
  expect_true(any(r1$history$event == "wgd"))

  # BFB detection and plots still work with WGD nodes (NA deltas)
  expect_s3_class(detect_bfb_branches(r1), "data.frame")
  expect_true(all(is.finite(r1$tree$edge.length)))
})

test_that("fit(wgd = TRUE) changes nothing when no cell is doubled", {
  sim <- wgd_sim(5, 0)
  chrs <- as.character(1:4)
  r0 <- fit(sim$cna_data, chromosomes = chrs, verbose = FALSE)
  r1 <- fit(sim$cna_data, chromosomes = chrs, wgd = TRUE, verbose = FALSE)
  expect_identical(r1$D, r0$D)
  expect_equal(nrow(r1$wgd$edges), 0)
  expect_equal(r1$reconstructions, r0$reconstructions)
})

test_that("restricting WGD costs to differing counts gives the same fit", {
  for (seed in c(2, 7)) {
    sim <- wgd_sim(seed, 2)
    chrs <- as.character(1:4)
    fast <- fit(sim$cna_data, chromosomes = chrs, wgd = TRUE, verbose = FALSE)
    full <- fit(sim$cna_data, chromosomes = chrs, wgd = TRUE, wgd_all_pairs = TRUE, verbose = FALSE)
    expect_identical(fast$D, full$D)
    expect_equal(fast$tree, full$tree)
    expect_equal(fast$wgd$edges$cells, full$wgd$edges$cells)
  }
})

test_that("fit(wgd = TRUE) validates its inputs", {
  sim <- wgd_sim(5, 0)
  expect_error(fit(sim$cna_data, chromosomes = "1", wgd = TRUE, g_dist_func = "greedy",
                   verbose = FALSE), "requires the default distances")
  expect_error(fit(sim$cna_data, chromosomes = "1", wgd = TRUE, wgd_counts = c(1, 2),
                   verbose = FALSE), "named vector")
  expect_error(fit(sim$cna_data, chromosomes = "1", wgd = TRUE, wgd_penalty = -1,
                   verbose = FALSE), "wgd_penalty")
})
