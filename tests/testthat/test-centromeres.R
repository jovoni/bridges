test_that("centromeres dataset has both builds and all 24 standard chromosomes", {
  expect_setequal(unique(centromeres$genome), c("hg19", "hg38"))

  main_chrs <- c(as.character(1:22), "X", "Y")
  expect_setequal(centromeres$chr[centromeres$genome == "hg19"], main_chrs)
  expect_setequal(centromeres$chr[centromeres$genome == "hg38"], main_chrs)
  expect_equal(nrow(centromeres), 2L * length(main_chrs))
})

test_that("centromeres coordinates are internally consistent", {
  expect_true(all(centromeres$centromere_start < centromeres$centromere_end))
  expect_true(all(centromeres$centromere_end < centromeres$chr_length))
  expect_true(all(centromeres$p_arm_start == 1L))
  expect_true(all(centromeres$p_arm_end == centromeres$centromere_start - 1L))
  expect_true(all(centromeres$q_arm_start == centromeres$centromere_end + 1L))
  expect_true(all(centromeres$q_arm_end == centromeres$chr_length))
  expect_false(anyNA(centromeres))
})

test_that("centromeres chr1 values match known real coordinates for hg19/hg38", {
  hg19_chr1 <- centromeres[centromeres$genome == "hg19" & centromeres$chr == "1", ]
  hg38_chr1 <- centromeres[centromeres$genome == "hg38" & centromeres$chr == "1", ]

  # UCSC-published values: hg19 chr1 length 249,250,621 bp, centromere
  # ~121.5-124.5 Mb; hg38 chr1 length 248,956,422 bp, centromere ~122.0-124.9 Mb.
  expect_equal(hg19_chr1$chr_length, 249250621L)
  expect_equal(hg19_chr1$centromere_start, 121535435L, tolerance = 1e4)
  expect_equal(hg38_chr1$chr_length, 248956422L)
  expect_equal(hg38_chr1$centromere_start, 122026460L, tolerance = 1e4)
})

test_that(".chr_lengths_for resolves real lengths and errors on unknown chromosomes", {
  lens <- bridges:::.chr_lengths_for("hg19", c("1", "2"))
  expect_equal(unname(lens), c(249250621, 243199373))

  expect_error(bridges:::.chr_lengths_for("hg19", "MT"),
               "No chr_length available")
  expect_error(bridges:::.chr_lengths_for("hg37", "1"), "should be one of")
})
