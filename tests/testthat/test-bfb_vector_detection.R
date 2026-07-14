test_that("bfb_detect_vectors: correctness + dedup", {
  skip_if_not_installed("bfbtools")
  vec_list <- list(
    c(6, 3, 5), c(6, 3, 5),               # exact duplicate pair
    c(14, 7, 18, 16, 9, 12),              # non-admitting, near-miss example
    NULL,                                  # uninformative
    c(5),                                  # length-1, uninformative
    c(NA, 5, 3),                           # NA-containing, uninformative
    c(14, 7, 18, 16, 9, 12),               # duplicate of row 3
    c(1, 1),                               # trivial: flat 1s + deletion gap
    c(1, 1, 1, 1)                          # trivial: same idea, longer
  )
  res <- bfb_detect_vectors(vec_list)

  expect_true(res$admits_bfb[1]); expect_equal(res$bfb_distance[1], 0)
  expect_true(res$admits_bfb[2])
  expect_false(res$admits_bfb[3])
  expect_equal(res$bfb_distance[3], -log(0.9114639), tolerance = 1e-4)
  expect_true(is.na(res$admits_bfb[4]))
  expect_true(is.na(res$admits_bfb[5]))
  expect_true(is.na(res$admits_bfb[6]))
  expect_false(res$admits_bfb[7])
  expect_true(is.na(res$admits_bfb[8]))  # trivial vector correctly excluded, not TRUE
  expect_true(is.na(res$admits_bfb[9]))  # same
  expect_equal(res$n_unique[1], 2)       # only 2 distinct informative vectors among 9 rows
})

test_that("split_chr_arms labels p/q correctly", {
  arm_test_data <- rbind(
    data.frame(chr = 1, start = c(0, 10, 20, 30, 40, 50), end = c(9,19,29,39,49,59)),
    data.frame(chr = 2, start = c(0, 100), end = c(99, 199))
  )
  arm_out <- split_chr_arms(arm_test_data, centromere = c("1" = 30, "2" = 50))

  expect_true(all(arm_out$arm[arm_out$chr == 1 & arm_out$start < 30] == "p"))
  expect_true(all(arm_out$arm[arm_out$chr == 1 & arm_out$start >= 30] == "q"))
  expect_true(all(arm_out$arm[arm_out$chr == 2 & arm_out$start < 50] == "p"))
  expect_true(all(arm_out$arm[arm_out$chr == 2 & arm_out$start >= 50] == "q"))
})

test_that("bfb_detect_batch: whole-chromosome mode end to end", {
  skip_if_not_installed("bfbtools")
  make_rows <- function(cell_id, chr, A, B) {
    n <- length(A)
    data.frame(cell_id = cell_id, chr = chr,
               start = seq(0, by = 10, length.out = n),
               end = seq(9, by = 10, length.out = n),
               A = A, B = B)
  }
  cna_data <- rbind(
    make_rows("cellX", 1, rev(c(6, 3, 5)), c(1, 1, 1)),
    make_rows("cellX", 2, rev(c(14, 7, 18, 16, 9, 12)), c(2, 2, 2)),
    make_rows("cellY", 1, rev(c(6, 3, 5)), c(1, 1, 1)),
    make_rows("cellY", 2, c(1, 1), c(1, 1))
  )
  out <- bfb_detect_batch(cna_data, alleles = c("A", "B"))

  expect_setequal(unique(out$cell_id), c("cellX", "cellY"))
  expect_equal(nrow(out), 8)
  expect_true(out$admits_bfb[out$cell_id == "cellX" & out$chr == 1 & out$allele == "A"])
  expect_true(is.na(out$admits_bfb[out$cell_id == "cellX" & out$chr == 1 & out$allele == "B"]))
  expect_true(out$admits_bfb[out$cell_id == "cellY" & out$chr == 1 & out$allele == "A"])
  expect_false(out$admits_bfb[out$cell_id == "cellX" & out$chr == 2 & out$allele == "A"])
})

test_that("bfb_detect_batch: arm-aware mode recovers the true pattern on both arms", {
  skip_if_not_installed("bfbtools")
  # cellZ carries the true [6,3,5] BFB pattern on BOTH arms of chr 1, laid
  # out in raw ascending-coordinate order per arm's own convention: p-arm
  # ascending = telomere-first already (raw = true vector directly);
  # q-arm ascending = centromere-first (raw = reversed true vector). A
  # whole-chromosome (non-arm-aware) reading of this row would see a
  # single 6-segment run and get it wrong.
  arm_cna <- rbind(
    data.frame(cell_id = "cellZ", chr = 1, start = c(0,10,20), end = c(9,19,29),
               A = c(6,3,5), B = c(1,1,1)),        # p-arm: raw = true vector
    data.frame(cell_id = "cellZ", chr = 1, start = c(30,40,50), end = c(39,49,59),
               A = c(5,3,6), B = c(1,1,1))         # q-arm: raw = reversed true vector
  )
  arm_result <- bfb_detect_batch(arm_cna, alleles = "A", centromere = c("1" = 30))

  expect_equal(nrow(arm_result), 2)  # p and q, allele A only
  expect_true(all(arm_result$admits_bfb))
})
