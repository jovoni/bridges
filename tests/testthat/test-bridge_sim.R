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

test_that("bridge_sim(genome=) wires through to real chr_lengths and bin counts", {
  set.seed(5)
  sim_hg19 <- bridge_sim(initial_cells = 1, chromosomes = c("1"), genome = "hg19",
                          bin_length = 1e4,
                          bfb_prob = 0.2, amp_rate = 0.2, del_rate = 0.05,
                          max_cells = 10, lambda = 1)
  set.seed(5)
  sim_hg38 <- bridge_sim(initial_cells = 1, chromosomes = c("1"), genome = "hg38",
                          bin_length = 1e4,
                          bfb_prob = 0.2, amp_rate = 0.2, del_rate = 0.05,
                          max_cells = 10, lambda = 1)

  expect_equal(sim_hg19$input_parameters$genome, "hg19")
  expect_equal(unname(sim_hg19$input_parameters$chr_lengths["1"]), 249250621)
  expect_equal(unname(sim_hg38$input_parameters$chr_lengths["1"]), 248956422)
  expect_false(identical(sim_hg19$input_parameters$chr_seq_lengths,
                          sim_hg38$input_parameters$chr_seq_lengths))
})

test_that("bridge_sim errors clearly on an invalid genome", {
  expect_error(
    bridge_sim(initial_cells = 1, chromosomes = c("1"), genome = "hg37", max_cells = 5),
    "should be one of"
  )
})

test_that("bridge_sim supports custom breakpoints throughout the simulation", {
  set.seed(3)
  sim <- bridge_sim(initial_cells = 1, chromosomes = c("1"), bfb_prob = 1,
                     amp_rate = 0, del_rate = 0, max_cells = 20, lambda = 2,
                     breakpoint_support = "custom", custom_breakpoints = c(50, 150, 200))
  expect_true(is.data.frame(sim$cna_data))
  expect_gt(sum(sim$cell_history$cn_event != "none"), 2)

  # On an unaltered chromosome (1..249, fused length 498) the cut must fall
  # right after a custom bin v: at position v on the forward half (shorter
  # daughter v bins) or at 499 - v on the reversed half (shorter daughter v - 1).
  seq0 <- bridges:::vec2seq(1:249)
  for (i in 1:20) {
    d <- bridges:::sim_bfb_left_and_right_sequences(
      seq0, support = "custom", custom_breakpoints = c(50, 150, 200))
    lens <- sort(c(bridges:::get_seq_length(d$l_seq), bridges:::get_seq_length(d$r_seq)))
    expect_true(lens[1] %in% c(50, 150, 200, 49, 149, 199))
  }
})

test_that("bridge_sim seeded from initial_sequences returns a valid tree", {
  set.seed(7)
  s0 <- bridge_sim(chromosomes = "1", max_cells = 10, return_phylo = FALSE)
  expect_no_warning(
    s1 <- bridge_sim(initial_sequences = s0$cells, chromosomes = "1", max_cells = 20)
  )
  expect_setequal(s1$tree$tip.label, names(s1$cells))
})

test_that("bridge_sim never reports an elapsed time beyond max_time", {
  set.seed(9)
  s <- bridge_sim(chromosomes = "1", max_cells = 1e4, max_time = 5, birth_rate = 1,
                  return_phylo = FALSE, return_cna_data = FALSE)
  expect_lte(s$elapsed_time, 5)
})

test_that("bridge_sim validates lambda, rate, subsample and hotspot position", {
  expect_error(bridge_sim(chromosomes = "1", lambda = -1), "lambda")
  expect_error(bridge_sim(chromosomes = "1", rate = 0), "rate")
  expect_error(bridge_sim(chromosomes = "1", subsample = 0), "subsample")
  expect_error(bridge_sim(chromosomes = "1", breakpoint_support = "custom"), "custom_breakpoints")
  expect_warning(
    bridge_sim(chromosomes = "21", bfb_allele = "21:A",
               hotspot = list(chr = "21:A", pos = 500), max_cells = 4),
    "beyond the last bin"
  )
})
