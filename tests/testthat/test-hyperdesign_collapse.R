library(testthat)
library(tibble)

test_that("as_multidesign.hyperdesign collapses correctly", {
  d1 <- multidesign(
    matrix(1:20, 10, 2),
    data.frame(condition = rep(c("A", "B"), 5))
  )
  d2 <- multidesign(
    matrix(21:40, 10, 2),
    data.frame(condition = rep(c("A", "B"), 5))
  )
  hd <- hyperdesign(list(d1, d2))

  md <- as_multidesign(hd)
  expect_s3_class(md, "multidesign")
  expect_equal(nrow(md$x), 20)
  expect_equal(ncol(md$x), 2)
  expect_equal(nrow(md$design), 20)
  expect_true("condition" %in% names(design(md)))
})

test_that("as_multidesign.hyperdesign adds .id column", {
  d1 <- multidesign(
    matrix(rnorm(10), 5, 2),
    data.frame(cond = rep("A", 5))
  )
  d2 <- multidesign(
    matrix(rnorm(10), 5, 2),
    data.frame(cond = rep("B", 5))
  )
  hd <- hyperdesign(list(d1, d2), block_names = c("block1", "block2"))

  md <- as_multidesign(hd, .id = "source")
  expect_true("source" %in% names(design(md)))
  expect_equal(md$design$source, c(rep("block1", 5), rep("block2", 5)))
})

test_that("as_multidesign.hyperdesign errors on mismatched ncol", {
  d1 <- multidesign(
    matrix(rnorm(10), 5, 2),
    data.frame(cond = rep("A", 5))
  )
  d2 <- multidesign(
    matrix(rnorm(15), 5, 3),
    data.frame(cond = rep("B", 5))
  )
  hd <- hyperdesign(list(d1, d2))

  expect_error(
    as_multidesign(hd),
    "all blocks must have the same number of columns"
  )
})

test_that("as_multidesign.hyperdesign errors on mismatched column_design", {
  d1 <- multidesign(
    matrix(rnorm(10), 5, 2),
    data.frame(cond = rep("A", 5)),
    data.frame(var = c("x", "y"))
  )
  d2 <- multidesign(
    matrix(rnorm(10), 5, 2),
    data.frame(cond = rep("B", 5)),
    data.frame(var = c("a", "b"))
  )
  hd <- hyperdesign(list(d1, d2))

  expect_error(
    as_multidesign(hd),
    "column designs must be identical"
  )
})

test_that("bind_multidesign handles mix of hyperdesign and multidesign", {
  col_design <- tibble(feature = 1:3)
  d1 <- multidesign(
    matrix(rnorm(15), 5, 3),
    data.frame(cond = rep("A", 5)),
    col_design
  )
  d2 <- multidesign(
    matrix(rnorm(15), 5, 3),
    data.frame(cond = rep("B", 5)),
    col_design
  )
  hd <- hyperdesign(list(d1, d2))

  d3 <- multidesign(
    matrix(rnorm(15), 5, 3),
    data.frame(cond = rep("C", 5)),
    col_design
  )

  combined <- bind_multidesign(hd, d3)
  expect_s3_class(combined, "multidesign")
  expect_equal(nrow(combined$x), 15)
  expect_equal(ncol(combined$x), 3)
  expect_equal(nrow(combined$design), 15)
})
