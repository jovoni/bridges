# Checks for the event-model kernels (dev/wgd/event_model.cpp).
#   Rscript dev/wgd/check_event_model.R
Rcpp::sourceCpp("dev/wgd/event_model.cpp")
ok <- function(label, cond) cat(sprintf("%-70s %s\n", label, if (isTRUE(cond)) "OK" else "FAILED"))

# 1. E(a, b) equals the true minimum number of +-1 segment events (BFS) ----
min_events_bfs <- function(a, b, lo = -3, hi = 7) {
  n <- length(a); key <- function(v) paste(v, collapse = ",")
  seen <- setNames(0L, key(a)); frontier <- list(a); steps <- 0L
  while (length(frontier)) {
    for (v in frontier) if (all(v == b)) return(steps)
    steps <- steps + 1L; nxt <- list()
    for (v in frontier) for (i in 1:n) for (j in i:n) for (s in c(-1, 1)) {
      w <- v; w[i:j] <- w[i:j] + s
      if (any(w < lo | w > hi)) next
      k <- key(w); if (!is.null(seen[k]) && !is.na(seen[k])) next
      seen[k] <- steps; nxt[[length(nxt) + 1]] <- w
    }
    frontier <- nxt
  }
  NA
}
set.seed(1)
bfs_ok <- replicate(150, {
  n <- sample(1:4, 1); a <- sample(0:3, n, TRUE); b <- sample(0:3, n, TRUE)
  event_distance(a, b) == min_events_bfs(a, b)
})
ok(sprintf("E(a, b) == brute-force minimum events (%d random pairs, n <= 4)", length(bfs_ok)), all(bfs_ok))

# 2. Ancestor DP is exact: compare with exhaustive search over x ------------
E <- function(a, b) event_distance(a, b)
exact_ok <- replicate(300, {
  n <- sample(1:4, 1); f <- sample(list(c(1, 1), c(1, 2), c(1, 4), c(2, 2), c(2, 1)), 1)[[1]]
  p <- sample(0:6, n, TRUE); q <- sample(0:9, n, TRUE)
  grid <- as.matrix(expand.grid(rep(list(0:10), n)))
  brute <- min(apply(grid, 1, function(x) E(f[1] * x, p) + E(f[2] * x, q)))
  r <- event_ancestor(p, f[1], q, f[2])
  isTRUE(all.equal(r$cost, brute)) && isTRUE(all.equal(r$cost, r$cost_p + r$cost_q))
})
ok(sprintf("ancestor DP == exhaustive search over x (%d random cases)", length(exact_ok)), all(exact_ok))

# 3. Same-level ancestor lies on a shortest path ----------------------------
path_ok <- replicate(300, {
  n <- sample(2:30, 1); a <- sample(0:5, n, TRUE); b <- sample(0:5, n, TRUE)
  r <- event_ancestor(a, 1, b, 1, t = 1)
  r$cost == E(a, b)
})
ok("same-level ancestor: E(x, a) + E(x, b) == E(a, b)", all(path_ok))

# 4. Known cases -------------------------------------------------------------
x <- c(1, 1, 3, 3, 1, 0, 0, 2, 2, 1)
ok("x -> 2x costs 0 events (+ w)", event_ancestor(x, 1, 2 * x, 2)$cost == 0)
ok("x -> 2x then one +1 gain costs 1", event_ancestor(x, 1, replace(2 * x, 3:4, 7), 2)$cost == 1)
ok("x -> 4x costs 0", event_ancestor(x, 1, 4 * x, 4)$cost == 0)
ok("2 -> 8 on one segment costs 6 events (greedy G counts 1)", event_distance(c(2, 2, 2), c(2, 8, 2)) == 6)
cat(sprintf("for reference: E(x, 2x) without WGD = %g events\n", E(x, 2 * x)))
