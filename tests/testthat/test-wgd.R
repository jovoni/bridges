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
  m <- bridges:::wgd_matrices_cpp(X)
  for (i in 1:20) for (j in 1:20) {
    if (i == j) next
    expect_equal(m$W1[i, j], W(X[i, ], X[j, ], 2))
    expect_equal(m$W2[i, j], W(X[i, ], X[j, ], 4))
  }
  expect_true(all(is.infinite(diag(m$W1))))
  expect_equal(m$n_unique, nrow(unique(X)))
})
