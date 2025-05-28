library(testthat)

# Test for orientation detection in multiblock

test_that("multiblock correctly detects orientation", {
  A1 <- matrix(1:6, nrow = 2)
  A2 <- matrix(7:12, nrow = 2)
  mb_c <- multiblock(list(A1, A2))
  expect_true(is_cstacked(mb_c))
  expect_false(is_rstacked(mb_c))

  B1 <- matrix(1:6, ncol = 2)
  B2 <- matrix(7:12, ncol = 2)
  mb_r <- multiblock(list(B1, B2))
  expect_true(is_rstacked(mb_r))
  expect_false(is_cstacked(mb_r))
})
