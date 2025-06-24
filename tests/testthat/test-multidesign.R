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