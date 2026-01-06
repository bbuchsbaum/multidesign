test_that("multidesign object creation works correctly", {
  # Create test data
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), c(2, 3)),
    subject = 1:5
  )
  col_info <- tibble(
    region = paste0("roi_", 1:4),
    type = rep(c("cortical", "subcortical"), 2)
  )
  
  # Test basic creation
  md <- multidesign(X, Y)
  expect_s3_class(md, "multidesign")
  expect_equal(dim(md$x), c(5, 4))
  expect_equal(nrow(md$design), 5)
  expect_true(!is.null(md$column_design))
  
  # Test with column design
  md_col <- multidesign(X, Y, col_info)
  expect_equal(md_col$column_design, col_info)
  
  # Test error conditions
  expect_error(
    multidesign(X[1:3,], Y),
    "`nrow\\(x\\)` must be equal to 5L"
  )
  expect_error(
    multidesign(X, Y, col_info[1:2,]),
    "`ncol\\(x\\)` must be equal to 2L"
  )
})

test_that("subsetting works correctly", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), c(2, 3)),
    subject = 1:5
  )
  md <- multidesign(X, Y)
  
  # Test simple subsetting
  sub_A <- subset(md, condition == "A")
  expect_equal(nrow(sub_A$x), 2)
  expect_equal(sub_A$design$condition, c("A", "A"))
  
  # Test multiple conditions
  sub_AB1 <- subset(md, condition == "A" & subject <= 2)
  expect_equal(nrow(sub_AB1$x), 2)
  
  # Test empty result
  sub_none <- subset(md, condition == "C")
  expect_null(sub_none)
})

test_that("splitting works correctly", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), c(2, 3)),
    subject = rep(1:3, length.out=5)
  )
  md <- multidesign(X, Y)
  
  # Test splitting by one variable
  split_cond <- split(md, condition)
  expect_length(split_cond, 2)  # A and B groups
  expect_equal(nrow(split_cond[[1]]$x), 2)  # A group
  expect_equal(nrow(split_cond[[2]]$x), 3)  # B group
  
  # Test splitting by multiple variables
  split_both <- split(md, condition, subject)
  expect_true(length(split_both) >= 2)
  
  # Test split_indices
  indices <- split_indices(md, condition)
  expect_equal(length(indices$indices), 2)
})

test_that("summarization works correctly", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), c(2, 3)),
    subject = rep(1:3, length.out=5)
  )
  md <- multidesign(X, Y)
  
  # Test summarization by condition
  sum_cond <- summarize_by(md, condition)
  expect_equal(nrow(sum_cond$x), 2)  # One row per condition
  expect_equal(ncol(sum_cond$x), ncol(X))
  expect_equal(sum_cond$design$condition, c("A", "B"))
  
  # Test custom summary function
  sum_custom <- summarize_by(md, condition, sfun=function(x) apply(x, 2, sd))
  expect_equal(dim(sum_custom$x), dim(sum_cond$x))
})

test_that("data extraction methods work", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), c(2, 3)),
    subject = 1:5
  )
  md <- multidesign(X, Y)

  # Test xdata extraction
  expect_equal(xdata(md), X)

  # Test design extraction
  expect_equal(design(md), Y)

  # Test column_design extraction
  expect_true(!is.null(column_design(md)))
})

# ============================================================================
# Regression tests for bug fixes
# ============================================================================

test_that("split.multidesign returns named list with correct names", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), c(2, 3)),
    subject = 1:5
  )
  md <- multidesign(X, Y)

  # Test that split returns a named list
  split_cond <- split(md, condition)
  expect_true(!is.null(names(split_cond)))
  expect_true("A" %in% names(split_cond))
  expect_true("B" %in% names(split_cond))

  # Test multi-variable split names (should be underscore-separated)
  X2 <- matrix(1:24, 6, 4)
  Y2 <- tibble(
    cond = rep(c("X", "Y"), each = 3),
    block = rep(1:3, 2)
  )
  md2 <- multidesign(X2, Y2)
  split_both <- split(md2, cond, block)
  expect_true(!is.null(names(split_both)))
  # Names should be like "X_1", "X_2", etc.
  expect_true(all(grepl("_", names(split_both))))
})

test_that("split.multidesign preserves correct data using .row_id indexing", {
  # Create data where row order matters
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = c("B", "A", "B", "A", "B"),  # Non-contiguous groups
    value = c(10, 20, 30, 40, 50)
  )
  md <- multidesign(X, Y)

  split_result <- split(md, condition)

  # Verify A group has correct rows (rows 2 and 4)
  A_group <- split_result[["A"]]
  expect_equal(nrow(A_group$x), 2)
  expect_equal(A_group$design$value, c(20, 40))
  # Check that the actual data rows are correct
  expect_equal(A_group$x[1, 1], X[2, 1])  # Row 2 of original
  expect_equal(A_group$x[2, 1], X[4, 1])  # Row 4 of original

  # Verify B group has correct rows (rows 1, 3, 5)
  B_group <- split_result[["B"]]
  expect_equal(nrow(B_group$x), 3)
  expect_equal(B_group$design$value, c(10, 30, 50))
  expect_equal(B_group$x[1, 1], X[1, 1])  # Row 1 of original
  expect_equal(B_group$x[2, 1], X[3, 1])  # Row 3 of original
  expect_equal(B_group$x[3, 1], X[5, 1])  # Row 5 of original
})

test_that("split.multidesign preserves column_design", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(condition = rep(c("A", "B"), c(2, 3)))
  col_info <- tibble(
    region = paste0("roi_", 1:4),
    type = rep(c("cortical", "subcortical"), 2)
  )
  md <- multidesign(X, Y, col_info)

  split_result <- split(md, condition)

  # Both splits should have the same column_design
  expect_equal(split_result[["A"]]$column_design, col_info)
  expect_equal(split_result[["B"]]$column_design, col_info)
})

test_that("split_indices.multidesign uses correct .row_id indexing", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = c("B", "A", "B", "A", "B"),  # Non-contiguous groups
    value = 1:5
  )
  md <- multidesign(X, Y)

  indices <- split_indices(md, condition)

  # A should have indices 2 and 4
  A_indices <- indices$indices[indices$condition == "A"][[1]]
  expect_equal(sort(A_indices), c(2, 4))

  # B should have indices 1, 3, 5
  B_indices <- indices$indices[indices$condition == "B"][[1]]
  expect_equal(sort(B_indices), c(1, 3, 5))
})

test_that("summarize_by.multidesign uses correct .row_id indexing", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = c("B", "A", "B", "A", "B"),  # Non-contiguous groups
    value = 1:5
  )
  md <- multidesign(X, Y)

  # Calculate expected means manually
  # A group: rows 2 and 4
  expected_A_means <- colMeans(X[c(2, 4), , drop = FALSE])
  # B group: rows 1, 3, 5
  expected_B_means <- colMeans(X[c(1, 3, 5), , drop = FALSE])

  result <- summarize_by(md, condition)

  # Find which row is A and which is B
  A_row <- which(result$design$condition == "A")
  B_row <- which(result$design$condition == "B")

  expect_equal(as.vector(result$x[A_row, ]), expected_A_means)
  expect_equal(as.vector(result$x[B_row, ]), expected_B_means)
})

test_that("summarize_by.multidesign preserves column_design", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(condition = rep(c("A", "B"), c(2, 3)))
  col_info <- tibble(
    region = paste0("roi_", 1:4),
    type = rep(c("cortical", "subcortical"), 2)
  )
  md <- multidesign(X, Y, col_info)

  result <- summarize_by(md, condition)
  expect_equal(result$column_design, col_info)
})

test_that("fold_over.multidesign preserves column_design", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(condition = rep(c("A", "B"), c(2, 3)))
  col_info <- tibble(
    region = paste0("roi_", 1:4),
    type = rep(c("cortical", "subcortical"), 2)
  )
  md <- multidesign(X, Y, col_info)

  folds <- fold_over(md, condition)

  # Check first fold
  fold1 <- folds[[1]]
  expect_equal(fold1$analysis$column_design, col_info)
  expect_equal(fold1$assessment$column_design, col_info)
})

test_that("fold_over.multidesign handles single-row folds correctly with drop=FALSE", {
  X <- matrix(1:12, 3, 4)
  Y <- tibble(condition = c("A", "B", "C"))
  md <- multidesign(X, Y)

  folds <- fold_over(md, condition)

  # Each assessment set should have 1 row but still be a matrix
  for (i in seq_along(folds)) {
    fold <- folds[[i]]
    expect_true(is.matrix(fold$assessment$x))
    expect_equal(nrow(fold$assessment$x), 1)
    expect_equal(ncol(fold$assessment$x), 4)
  }
})

test_that("reduce generic exists and reduce.multidesign is S3 method",
{
  # Test that reduce generic exists
  expect_true(exists("reduce"))
  expect_true(is.function(reduce))

  # Test that reduce.multidesign works via S3 dispatch
  X <- matrix(rnorm(50), 10, 5)
  Y <- tibble(condition = rep(c("A", "B"), each = 5))
  md <- multidesign(X, Y)

  # This should work via S3 dispatch
  reduced <- reduce(md, nc = 2)
  expect_s3_class(reduced, "reduced_multidesign")
  expect_s3_class(reduced, "multidesign")
  expect_equal(ncol(reduced$x), 2)
})

test_that("print.reduced_multidesign handles missing column_design", {
  X <- matrix(rnorm(50), 10, 5)
  Y <- tibble(condition = rep(c("A", "B"), each = 5))
  md <- multidesign(X, Y)

  reduced <- reduce(md, nc = 2)

  # reduced_multidesign does not have column_design
  expect_null(reduced$column_design)

 # print should work without error
  expect_output(print(reduced), "Reduced Multidesign Object")
})

test_that("print.multidesign works correctly", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), c(2, 3)),
    subject = 1:5
  )
  md <- multidesign(X, Y)

  # Should print without error
  expect_output(print(md), "Multidesign Object")
  expect_output(print(md), "5 observations")
  expect_output(print(md), "4 variables")
})