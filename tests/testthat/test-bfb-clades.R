# A minimal fit-like object: a balanced rooted tree and hand-set deltas for
# one chromosome/allele (positive = BFB-like node)
toy_fit <- function(n_tip, bfb_nodes) {
  tree <- ape::stree(n_tip, "balanced")
  internal <- n_tip + seq_len(tree$Nnode)
  deltas <- stats::setNames(as.list(ifelse(internal %in% bfb_nodes, 1, -1)), internal)
  list(tree = tree,
       all_input_Xs = list("1" = list(A = matrix(0, n_tip, 1, dimnames = list(tree$tip.label, NULL)))),
       reconstructions = list("1" = list(A = list(deltas = deltas))))
}
descendant_internal <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  sub <- ape::extract.clade(tree, node)
  node:(node + sub$Nnode - 1)   # ape numbers the internal nodes of a clade contiguously
}

test_that("BFB-like nodes packed in one clade are found even when diluted genome-wide", {
  f <- toy_fit(512, bfb_nodes = NULL)
  # a clade of 15 internal nodes (16 tips)
  clade <- which(sapply(seq_len(f$tree$Nnode) + 512, function(n) length(ape::extract.clade(f$tree, n)$tip.label)) == 16)[1] + 512
  packed <- descendant_internal(f$tree, clade)[c(1, 3, 6, 9)]
  f <- toy_fit(512, packed)

  global <- detect_bfb_branches(f)
  expect_gt(global$p.value, 0.05)   # 4 of 511 nodes: diluted

  res <- detect_bfb_clades(f, n_sim = 499, seed = 1)
  expect_lt(res$p.scan, 0.01)
  expect_equal(res$total_bfb_nodes, 4)
  expect_equal(res$clade_bfb_nodes, 4)
  expect_true(all(res$clade_cells[[1]] %in% ape::extract.clade(f$tree, clade)$tip.label))
})

test_that("the same number of BFB-like nodes scattered over the tree is not significant", {
  f <- toy_fit(512, bfb_nodes = NULL)
  sizes <- sapply(seq_len(f$tree$Nnode) + 512, function(n) length(ape::extract.clade(f$tree, n)$tip.label))
  deep <- which(sizes == 128) + 512                    # four disjoint quarters of the tree
  scattered <- sapply(deep, function(n) descendant_internal(f$tree, n)[60])
  res <- detect_bfb_clades(toy_fit(512, scattered), n_sim = 499, seed = 1)
  expect_gt(res$p.scan, 0.05)
})

test_that("no BFB-like node gives p.scan = 1, and results are reproducible", {
  expect_equal(detect_bfb_clades(toy_fit(64, NULL), n_sim = 99, seed = 2)$p.scan, 1)

  f <- toy_fit(128, bfb_nodes = c(130, 131, 132))
  set.seed(10); before <- runif(1); set.seed(10)
  a <- detect_bfb_clades(f, n_sim = 199, seed = 3)
  after <- runif(1)
  b <- detect_bfb_clades(f, n_sim = 199, seed = 3)
  expect_equal(a, b)
  expect_equal(before, after)   # the session's random-number state is restored
})

test_that("detect_bfb_clades runs on a fit with WGD and skips WGD nodes", {
  set.seed(4)
  sim <- suppressWarnings(bridge_sim(chromosomes = c("1", "2"), bfb_allele = "1:A", lambda = 1,
                                     max_cells = 120, subsample = 60, wgd_available = 1,
                                     wgd_probability = 0.1))
  res <- fit(sim$cna_data, chromosomes = c("1", "2"), wgd = TRUE, verbose = FALSE)
  out <- detect_bfb_clades(res, n_sim = 99, seed = 1)
  expect_equal(nrow(out), 4)
  expect_true(all(out$p.scan > 0 & out$p.scan <= 1))
  glob <- detect_bfb_branches(res)
  expect_equal(out$total_nodes, glob$N_trials[match(paste(out$chr, out$allele), paste(glob$chr, glob$allele))])
})

test_that("detect_bfb_clades validates its inputs", {
  f <- toy_fit(16, NULL)
  expect_error(detect_bfb_clades(f, threshold = 0), "threshold")
  expect_error(detect_bfb_clades(f, min_nodes = 0), "min_nodes")
  expect_error(detect_bfb_clades(f, n_sim = 0), "n_sim")
})
