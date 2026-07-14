test_that("bridge_sim runs and produces the expected structure", {
  set.seed(1)
  sim <- bridge_sim(initial_cells = 1, chromosomes = c("1"), bfb_prob = 0.3,
                     amp_rate = 0.3, del_rate = 0.1, max_cells = 20, lambda = 2)

  expect_type(sim, "list")
  expect_true(all(c("cells", "cell_history", "tree", "cna_data") %in% names(sim)))
  expect_true(is.data.frame(sim$cna_data))
  expect_true(all(c("cell_id", "chr", "A", "B", "CN", "start", "end") %in% names(sim$cna_data)))
  expect_true(inherits(sim$tree, "phylo"))
  expect_gt(nrow(sim$cna_data), 0)
})

test_that("bridge_sim respects max_cells as an upper bound", {
  set.seed(2)
  sim <- bridge_sim(initial_cells = 1, chromosomes = c("1"), bfb_prob = 0.2,
                     amp_rate = 0.2, del_rate = 0.05, max_cells = 15, lambda = 1)
  expect_lte(sim$n_alive, 15)
})

test_that("bridge_sim is reproducible under a fixed seed", {
  set.seed(42)
  sim1 <- bridge_sim(initial_cells = 1, chromosomes = c("1"), bfb_prob = 0.3,
                      amp_rate = 0.3, del_rate = 0.1, max_cells = 20, lambda = 2)
  set.seed(42)
  sim2 <- bridge_sim(initial_cells = 1, chromosomes = c("1"), bfb_prob = 0.3,
                      amp_rate = 0.3, del_rate = 0.1, max_cells = 20, lambda = 2)
  expect_identical(sim1$cna_data, sim2$cna_data)
})

test_that("bridge_sim works with custom breakpoint_support (exercises get_seq_length)", {
  # sim_bfb_left_and_right_sequences's "custom" branch calls
  # get_seq_length() directly rather than going through the C++ engine --
  # a real gap surfaced by R CMD check's static analysis during the
  # bfbgill extraction (get_seq_length was one of two duplicate
  # definitions in the original bridges package; only one was live).
  #
  # NOTE: the C++ main-loop engine explicitly rejects support="custom"
  # ("use support='uniform' or 'beta'") -- custom breakpoints are only
  # supported for the one-time initialization BFB event
  # (first_round_of_bfb = TRUE, the default), so bfb_prob must be 0 here
  # to keep the main loop from attempting another (unsupported) custom
  # BFB event later in the simulation.
  set.seed(3)
  sim <- bridge_sim(initial_cells = 1, chromosomes = c("1"), bfb_prob = 0,
                     amp_rate = 0.2, del_rate = 0.1, max_cells = 20, lambda = 2,
                     first_round_of_bfb = TRUE,
                     breakpoint_support = "custom", custom_breakpoints = c(50, 150, 250))
  expect_true(is.data.frame(sim$cna_data))
  expect_gt(nrow(sim$cna_data), 0)
})
