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
