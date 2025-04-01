
test_that("No error", {
  x = gillespie_sim(1, 100, death_rate = 0, bfb_prob = .05)
  cell_history_to_newick(x$cell_history)
})
