#
# test_that("test cell distances", {
#   set.seed(1234)
#   x1 = gillespie_sim(initial_cells = 1, death_rate = 0, max_cells = 100, bfb_prob = .05)
#   x2 = gillespie_sim(initial_cells = 1, death_rate = 0, max_cells = 100, bfb_prob = 0)
#
#   expect_equal(0, cells_distance(x1, x1, order = FALSE))
#   expect_true(cells_distance(x1, x2, order = FALSE) != 0)
#
#   x3 = gillespie_sim(initial_cells = 1, death_rate = 0, max_cells = 100, bfb_prob = 0, initial_sequence_length = 10)
#   expect_error(cells_distance(x1, x3, order = FALSE))
# })
