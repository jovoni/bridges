# Random copy-number vectors: mostly small values, some long equal runs, and
# some amplifications, to exercise multi-step blocks and parity changes.
random_cn_pair <- function(n) {
  base <- rep(sample(0:3, ceiling(n / 4), replace = TRUE), each = 4)[seq_len(n)]
  a <- base + rbinom(n, 1, 0.2) * sample(0:6, n, replace = TRUE)
  b <- base + rbinom(n, 1, 0.2) * sample(0:6, n, replace = TRUE)
  list(a = as.numeric(a), b = as.numeric(b))
}

test_that("greedy_distance_cpp matches the pure-R reference implementation", {
  set.seed(11)
  for (k in 1:2000) {
    v <- random_cn_pair(sample(1:40, 1))
    tv <- sample(c(1, 2), 1)
    fast <- bridges:::greedy_distance_cpp(v$a, v$b, target_val = tv)
    ref <- bridges:::find_greedy_distance_with_steps(v$a, v$b, target_val = tv)
    for (f in c("cost", "ancestor", "a_cost", "b_cost", "a_delta", "b_delta")) {
      expect_equal(fast[[f]], ref[[f]], info = paste(f, "case", k))
    }
  }
})

test_that("greedy_bfb_distance_cpp matches the pure-R reference implementation", {
  set.seed(12)
  for (k in 1:2000) {
    v <- random_cn_pair(sample(1:40, 1))
    pen <- sample(c(0, 1, 1.5), 1)
    fast <- bridges:::greedy_bfb_distance_cpp(v$a, v$b, penalty = pen)
    ref <- bridges:::find_greedy_distance_A_contig(v$a, v$b, penalty = pen)
    for (f in c("cost", "ancestor", "a_cost", "b_cost", "a_delta", "b_delta", "left", "right")) {
      expect_equal(fast[[f]], ref[[f]], info = paste(f, "case", k))
    }
  }
})

# Cells x segments matrix with many repeated profiles, as in real data
random_cn_matrix <- function(n_cells, n_seg, n_distinct) {
  profiles <- t(replicate(n_distinct, random_cn_pair(n_seg)$a))
  X <- profiles[sample(n_distinct, n_cells, replace = TRUE), , drop = FALSE]
  rownames(X) <- paste0("cell_", seq_len(n_cells))
  X
}

test_that("fast_distance_matrices matches the per-pair computation", {
  set.seed(13)
  for (k in 1:30) {
    X <- random_cn_matrix(n_cells = sample(2:40, 1), n_seg = sample(1:30, 1),
                          n_distinct = sample(1:15, 1))
    pen <- sample(c(0, 1.5, -1), 1)
    G <- bridges:::compute_distance_matrix(X, bridges:::greedy_distance_cpp, target_val = 1)
    B <- bridges:::compute_distance_matrix(X, bridges:::greedy_bfb_distance_cpp, penalty = pen)

    fast <- bridges:::fast_distance_matrices(X, use_bfb = TRUE, penalty = pen, keep_parts = TRUE)
    expect_identical(fast$greedy, G)
    expect_identical(fast$bfb, B)
    expect_identical(fast$min, pmin(G, B))
    expect_equal(fast$n_unique, nrow(unique(X)))

    greedy_only <- bridges:::fast_distance_matrices(X, use_bfb = FALSE)
    expect_identical(greedy_only$min, G)
    expect_null(greedy_only$greedy)
    expect_null(greedy_only$bfb)
  }
})

test_that("fast_distance_matrices handles all-identical rows and a single segment", {
  X <- matrix(c(2, 2, 2), ncol = 1, dimnames = list(c("a", "b", "c"), NULL))
  fast <- bridges:::fast_distance_matrices(X)
  expect_equal(fast$n_unique, 1L)
  expect_equal(unname(fast$min), matrix(c(Inf, 0, 0, 0, Inf, 0, 0, 0, Inf), 3))
})
