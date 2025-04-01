
test_that("hotspot test", {
  cell = vec2seq(c(1:100))
  expect_false(is_hotspot_gained(cell, NULL))

  cell = vec2seq(c(1:100, 100:85))
  expect_false(is_hotspot_gained(cell, 50))
  expect_true(is_hotspot_gained(cell, 90))
})
