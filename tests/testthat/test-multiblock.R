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

test_that("multiblock creation works correctly", {
  # Test column-stacked multiblock (same number of rows)
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(13:24, 4, 3)
  X3 <- matrix(25:36, 4, 3)
  mb_c <- multiblock(list(X1, X2, X3))
  
  expect_s3_class(mb_c, "multiblock")
  expect_s3_class(mb_c, "multiblock_c")
  expect_equal(attr(mb_c, "orient"), "cstacked")
  expect_equal(dim(attr(mb_c, "ind")), c(3, 2))  # 3 blocks, start/end indices
  
  # Test row-stacked multiblock (same number of columns)
  Y1 <- matrix(1:10, 2, 5)
  Y2 <- matrix(11:25, 3, 5)
  Y3 <- matrix(26:35, 2, 5)
  mb_r <- multiblock(list(Y1, Y2, Y3))
  
  expect_s3_class(mb_r, "multiblock")
  expect_s3_class(mb_r, "multiblock_r")
  expect_equal(attr(mb_r, "orient"), "rstacked")
  expect_equal(dim(attr(mb_r, "ind")), c(3, 2))  # 3 blocks, start/end indices
})

test_that("block_indices extraction works", {
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(13:24, 4, 3)
  X3 <- matrix(25:36, 4, 3)
  mb_c <- multiblock(list(X1, X2, X3))
  
  # Test indices for column-stacked blocks
  expect_equal(block_indices(mb_c, 1), 1:3)
  expect_equal(block_indices(mb_c, 2), 4:6)
  expect_equal(block_indices(mb_c, 3), 7:9)
  
  # Test indices for row-stacked blocks
  Y1 <- matrix(1:10, 2, 5)
  Y2 <- matrix(11:25, 3, 5)
  Y3 <- matrix(26:35, 2, 5)
  mb_r <- multiblock(list(Y1, Y2, Y3))
  
  expect_equal(block_indices(mb_r, 1), 1:2)
  expect_equal(block_indices(mb_r, 2), 3:5)
  expect_equal(block_indices(mb_r, 3), 6:7)
})

test_that("is_cstacked and is_rstacked work correctly", {
  # Column-stacked
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(13:24, 4, 3)
  mb_c <- multiblock(list(X1, X2))
  
  expect_true(is_cstacked(mb_c))
  expect_false(is_rstacked(mb_c))
  
  # Row-stacked
  Y1 <- matrix(1:10, 2, 5)
  Y2 <- matrix(11:25, 3, 5)
  mb_r <- multiblock(list(Y1, Y2))
  
  expect_false(is_cstacked(mb_r))
  expect_true(is_rstacked(mb_r))
})

test_that("transpose works correctly", {
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(13:24, 4, 3)
  mb_c <- multiblock(list(X1, X2))
  
  mb_t <- t(mb_c)
  expect_true(is_rstacked(mb_t))
  expect_false(is_cstacked(mb_t))
  expect_equal(nrow(mb_t[[1]]), 3)
  expect_equal(ncol(mb_t[[1]]), 4)
})

test_that("multiblock handles errors correctly", {
  # Empty list
  expect_error(multiblock(list()), "List is empty")
  
  # Non-matrix elements
  expect_error(
    multiblock(list(1:10, matrix(1:10, 2, 5))),
    "Not all elements are matrices"
  )
  
  # Mismatched dimensions
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(13:18, 2, 3)  # Different number of rows
  X3 <- matrix(19:24, 3, 2)  # Different number of columns
  expect_error(
    multiblock(list(X1, X2, X3)),
    "all matrices must share either row or column dimension"
  )
  
  # Invalid block index
  X1 <- matrix(1:12, 4, 3)
  X2 <- matrix(13:24, 4, 3)
  mb <- multiblock(list(X1, X2))
  expect_error(
    block_indices(mb, 3),
    "`i` must be between 1 and 2, not 3"
  )
})