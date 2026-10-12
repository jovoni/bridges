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
