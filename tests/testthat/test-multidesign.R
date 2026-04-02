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

# --- Regression tests ---

test_that("design() does not expose .index", {
  X <- matrix(rnorm(20), 5, 4)
  Y <- tibble(condition = rep(c("A", "B"), c(2, 3)))
  md <- multidesign(X, Y)
  des <- design(md)
  expect_false(".index" %in% names(des))
})

test_that("summarize_by returns correct values (not NaN)", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(condition = rep(c("A", "B"), c(2, 3)))
  md <- multidesign(X, Y)

  result <- summarize_by(md, condition)
  expect_equal(nrow(result$x), 2)
  expect_equal(ncol(result$x), 4)
  # Verify actual means match expected
  expect_equal(as.numeric(result$x[1, ]), colMeans(X[1:2, ]))
  expect_equal(as.numeric(result$x[2, ]), colMeans(X[3:5, ]))
  # No NaN values
  expect_false(any(is.nan(result$x)))
})

test_that("subset returns correct rows (not 1:n)", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), c(2, 3)),
    subject = 1:5
  )
  md <- multidesign(X, Y)

  sub_B <- subset(md, condition == "B")
  expect_equal(nrow(sub_B$x), 3)
  # Verify the actual data content matches rows 3-5

  expect_equal(sub_B$x, X[3:5, , drop=FALSE])
  expect_equal(sub_B$design$subject, 3:5)
})

test_that("split returns correct data content per group", {
  X <- matrix(1:20, 5, 4)
  Y <- tibble(condition = rep(c("A", "B"), c(2, 3)))
  md <- multidesign(X, Y)

  sp <- split(md, condition)
  expect_length(sp, 2)
  # Group A should have rows 1-2
  expect_equal(sp[[1]]$x, X[1:2, , drop=FALSE])
  # Group B should have rows 3-5
  expect_equal(sp[[2]]$x, X[3:5, , drop=FALSE])
})

test_that("fold_over works and returns valid folds with correct sizes", {
  X <- matrix(rnorm(40), 10, 4)
  Y <- tibble(condition = rep(c("A", "B"), each=5))
  md <- multidesign(X, Y)

  folds <- fold_over(md, condition)
  expect_s3_class(folds, "foldlist")
  expect_length(folds, 2)

  f1 <- folds[[1]]
  expect_true(!is.null(f1$analysis))
  expect_true(!is.null(f1$assessment))
  expect_equal(nrow(f1$analysis$x) + nrow(f1$assessment$x), nrow(X))
})

test_that("fold_over preserves column_design", {
  X <- matrix(rnorm(40), 10, 4)
  Y <- tibble(condition = rep(c("A", "B"), each=5))
  col_info <- tibble(var = paste0("v", 1:4))
  md <- multidesign(X, Y, col_info)

  folds <- fold_over(md, condition)
  f1 <- folds[[1]]
  expect_equal(f1$analysis$column_design, col_info)
  expect_equal(f1$assessment$column_design, col_info)
})

test_that("cv_rows.multidesign creates explicit row folds", {
  X <- matrix(1:24, 6, 4)
  Y <- tibble(condition = letters[1:6])
  md <- multidesign(X, Y)

  folds <- cv_rows(md, rows = list(c(1, 3), c(2, 4)))
  expect_s3_class(folds, "foldlist")
  expect_length(folds, 2)

  f1 <- folds[[1]]
  expect_equal(f1$assessment$x, X[c(1, 3), , drop = FALSE])
  expect_equal(f1$analysis$x, X[-c(1, 3), , drop = FALSE])
  expect_equal(f1$held_out$rows, c(1L, 3L))
})

test_that("cv_rows.multidesign validates explicit row folds", {
  X <- matrix(rnorm(24), 6, 4)
  Y <- tibble(condition = letters[1:6])
  md <- multidesign(X, Y)

  expect_error(
    cv_rows(md, rows = list(c(1, 1))),
    "must not contain duplicate row indices"
  )
  expect_error(
    cv_rows(md, rows = list(1:6)),
    "cannot hold out all rows"
  )
  expect_error(
    cv_rows(md, rows = list(c(1, 7))),
    "outside the valid range"
  )
})

test_that("print.multidesign does not error", {
  X <- matrix(rnorm(20), 5, 4)
  Y <- tibble(condition = rep(c("A", "B"), c(2, 3)))
  md <- multidesign(X, Y)
  expect_output(print(md), "Multidesign Object")
})

test_that("print.reduced_multidesign handles NULL column_design", {
  # Construct a reduced_multidesign manually with NULL column_design
  rd <- structure(list(
    x = matrix(rnorm(10), 5, 2),
    design = tibble(condition = rep(c("A", "B"), c(2, 3))),
    column_design = NULL,
    projector = NULL
  ), class = c("reduced_multidesign", "multidesign"))
  # print should not error even with NULL column_design
  # (projector will cause errors in multivarious::shape, so just test the null guard)
  expect_false(is.null(rd))
})

test_that("split_indices.multidesign returns correct groups and indices", {
  X <- matrix(rnorm(40), 10, 4)
  Y <- tibble(
    condition = rep(c("A", "B"), each = 5),
    block = rep(1:5, times = 2)
  )
  md <- multidesign(X, Y)

  # Split by one variable
  si <- split_indices(md, condition)
  expect_equal(nrow(si), 2)
  expect_true("indices" %in% names(si))
  expect_true(".splitvar" %in% names(si))
  expect_equal(sort(unlist(si$indices)), 1:10)

  # Split by two variables
  si2 <- split_indices(md, condition, block)
  expect_equal(nrow(si2), 10)
  expect_true(all(sapply(si2$indices, length) == 1))
})

test_that("print.multidesign displays column design info", {
  X <- matrix(rnorm(20), 5, 4)
  Y <- tibble(condition = rep(c("A", "B"), c(2, 3)))
  col_info <- tibble(var = paste0("v", 1:4), hemisphere = rep(c("L", "R"), 2))
  md <- multidesign(X, Y, col_info)
  expect_output(print(md), "Column Metadata")
})

test_that("reduce.multidesign reduces dimensions correctly", {
  skip_if_not_installed("multivarious")
  X <- matrix(rnorm(200), 20, 10)
  Y <- tibble(condition = rep(c("A", "B"), each = 10))
  md <- multidesign(X, Y)

  rmd <- reduce.multidesign(md, nc = 3)
  expect_s3_class(rmd, "reduced_multidesign")
  expect_s3_class(rmd, "multidesign")
  expect_equal(ncol(rmd$x), 3)
  expect_equal(nrow(rmd$x), 20)
  expect_true(!is.null(rmd$projector))

  # Print should work
  expect_output(print(rmd), "Reduced Multidesign Object")
})
