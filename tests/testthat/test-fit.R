test_that("fit runs without spurious warnings and returns a rooted tree", {
  set.seed(4)
  sim <- bridge_sim(chromosomes = "1", max_cells = 20)
  expect_no_warning(res <- fit(sim$cna_data, chromosomes = "1", verbose = FALSE))
  expect_s3_class(res$tree, "phylo")
  expect_setequal(res$tree$tip.label, unique(sim$cna_data$cell_id))
  expect_equal(unname(diag(res$D)), rep(0, nrow(res$D)))

  out <- detect_bfb_branches(res)
  expect_setequal(out$allele, c("A", "B"))
  expect_true(all(out$adj.pval >= 0 & out$adj.pval <= 1))
})

test_that("fit rejects unknown distance-function names up front", {
  set.seed(4)
  sim <- bridge_sim(chromosomes = "1", max_cells = 10)
  expect_error(fit(sim$cna_data, chromosomes = "1", g_dist_func = "nope"), "Unknown g_dist_func")
  expect_error(fit(sim$cna_data, chromosomes = "1", b_dist_func = "nope"), "Unknown b_dist_func")
})

test_that("tibble_to_matrix does not modify a data.table input by reference", {
  dt <- data.table::data.table(cell_id = c("a", "a", "b", "b"), start = c(1, 11, 1, 11),
                               end = c(10, 20, 10, 20), A = c(1, 2, 1, 1))
  cols <- names(dt)
  m <- bridges:::tibble_to_matrix(dt, value_column = "A")
  expect_identical(names(dt), cols)
  expect_equal(dim(m), c(2L, 2L))
})

fit_test_data <- function(seed, chrs = c("1", "2"), cells = 40) {
  set.seed(seed)
  suppressWarnings(bridge_sim(chromosomes = chrs, bfb_allele = paste0(chrs[1], ":A"),
                              hotspot = list(chr = paste0(chrs[1], ":A"), pos = 20),
                              max_cells = 2 * cells, subsample = cells, lambda = 1,
                              return_phylo = FALSE))$cna_data
}

test_that("fit's fast distance path gives the same result as the per-pair R path", {
  d <- fit_test_data(21)
  settings <- list(list(), list(bfb_penalty = 1.5), list(alleles = "CN"))
  for (st in settings) {
    fast <- do.call(fit, c(list(d, chromosomes = c("1", "2"), verbose = FALSE), st))
    slow <- do.call(fit, c(list(d, chromosomes = c("1", "2"), verbose = FALSE,
                                g_dist_func = "greedy", b_dist_func = "bfb"), st))
    expect_identical(fast$D, slow$D)
    expect_equal(fast$tree, slow$tree)
    expect_equal(lapply(fast$reconstructions, lapply, `[[`, "deltas"),
                 lapply(slow$reconstructions, lapply, `[[`, "deltas"))
  }

  # Greedy-only
  fast <- fit(d, chromosomes = c("1", "2"), b_dist_func = NULL, verbose = FALSE)
  slow <- fit(d, chromosomes = c("1", "2"), g_dist_func = "greedy", b_dist_func = NULL,
              verbose = FALSE)
  expect_identical(fast$D, slow$D)
  expect_equal(fast$tree, slow$tree)
})

test_that("fit keeps per-chromosome matrices only when keep_matrices = TRUE", {
  d <- fit_test_data(22)
  res <- fit(d, chromosomes = c("1", "2"), verbose = FALSE)
  expect_null(res$greedy_Ds)
  expect_null(res$avg_Ds)

  res <- fit(d, chromosomes = c("1", "2"), keep_matrices = TRUE, verbose = FALSE)
  X <- res$all_input_Xs
  expect_identical(res$greedy_Ds,
                   bridges:::compute_greedy_distances(X, c("1", "2"), c("A", "B"), "greedy_fast"))
  expect_identical(res$avg_Ds,
                   bridges:::compute_avg_distances(X, c("1", "2"), c("A", "B"), 0, "bfb_fast"))

  slow <- fit(d, chromosomes = c("1", "2"), g_dist_func = "greedy", b_dist_func = "bfb",
              keep_matrices = TRUE, verbose = FALSE)
  expect_equal(res$greedy_Ds, slow$greedy_Ds)
  expect_equal(res$avg_Ds, slow$avg_Ds)
})

test_that("compute_fast_distances aligns cells by name across chromosomes", {
  d <- fit_test_data(23)
  X <- bridges:::process_input_data(dplyr::bind_rows(d, bridges:::create_diploid_data(d)),
                                    c("1", "2"), c("A", "B"), 0, 0)
  ref <- bridges:::compute_fast_distances(X, c("1", "2"), c("A", "B"))$D

  shuffled <- X
  perm <- sample(nrow(X[["2"]][["A"]]))
  shuffled[["2"]] <- lapply(X[["2"]], function(m) m[perm, , drop = FALSE])
  expect_identical(bridges:::compute_fast_distances(shuffled, c("1", "2"), c("A", "B"))$D, ref)

  missing_cell <- X
  missing_cell[["2"]] <- lapply(X[["2"]], function(m) m[-1, , drop = FALSE])
  expect_error(bridges:::compute_fast_distances(missing_cell, c("1", "2"), c("A", "B")),
               "same cells")
})
