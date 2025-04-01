
test_that("Correct number of cells is returned", {
  n_cells = 123
  x = gillespie_sim(1, 100, max_cells = n_cells, birth_rate = 1, death_rate = 0)
  expect_equal(length(x$cells), n_cells)
})

test_that("All dead cells are correctly handled", {
  x = gillespie_sim(1, 100, death_rate = 1)
  expect_equal(length(x$cells), 0)
})
