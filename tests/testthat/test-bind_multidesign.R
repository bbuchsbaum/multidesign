library(testthat)
library(tibble)

# Test for successful binding of multidesign objects

test_that("bind_multidesign merges objects and adds id column", {
  X1 <- matrix(1:6, nrow = 3)
  X2 <- matrix(7:12, nrow = 3)
  design1 <- tibble(group = letters[1:3])
  design2 <- tibble(group = letters[4:6])
  col_design <- tibble(feature = 1:2)

  md1 <- multidesign(X1, design1, col_design)
  md2 <- multidesign(X2, design2, col_design)

  out <- bind_multidesign(md1, md2, .id = "src")

  expect_equal(nrow(out$x), 6)
  expect_equal(out$design$src, c(rep(1, 3), rep(2, 3)))
  expect_equal(out$column_design, col_design)
})

# Test for mismatched column designs

test_that("bind_multidesign fails with mismatched column designs", {
  X1 <- matrix(1:6, nrow = 3)
  X2 <- matrix(7:12, nrow = 3)
  design1 <- tibble(group = letters[1:3])
  design2 <- tibble(group = letters[4:6])
  col_design1 <- tibble(feature = 1:2)
  col_design2 <- tibble(feature = 1:3)

  md1 <- multidesign(X1, design1, col_design1)
  md2 <- multidesign(X2, design2, col_design2)

  expect_error(bind_multidesign(md1, md2),
               "column designs must be identical across multidesigns")
})
