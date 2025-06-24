test_that("observation_set creation works correctly", {
  # Test matrix input
  X_mat <- matrix(1:20, 5, 4)
  obs_mat <- obs_group(X_mat)
  expect_s3_class(obs_mat, "observation_set")
  expect_length(obs_mat, 5)
  expect_equal(as.vector(obs_mat[[1]]), as.vector(X_mat[1,]))  # Compare as vectors
  
  # Test list input
  X_list <- list(a=1:4, b=5:8, c=9:12)
  obs_list <- obs_group(X_list)
  expect_s3_class(obs_list, "observation_set")
  expect_length(obs_list, 3)
  expect_equal(obs_list[[1]], X_list[[1]])
  
  # Test custom indices
  ind <- c(10, 20, 30)
  obs_ind <- obs_group(X_list, ind=ind)
  expect_length(obs_ind, 3)
  
  # Test error conditions
  expect_error(obs_group(list()), "must not be empty")
  expect_error(obs_group(X_list, ind=1:4), "must be equal to")
})

test_that("multiframe list constructor works correctly", {
  # Create test data
  x <- list(
    matrix(1:12, 3, 4),
    matrix(13:24, 3, 4),
    matrix(25:36, 3, 4)
  )
  y <- data.frame(
    condition = c("A", "B", "C"),
    block = 1:3
  )
  
  # Test basic creation
  mf <- multiframe(x, y)
  expect_s3_class(mf, "multiframe")
  expect_equal(nrow(mf$design), 3)
  expect_true(all(c("condition", "block", ".obs", ".index") %in% names(mf$design)))
  expect_equal(length(mf$design$.obs), 3)  # Check number of observation functions
  
  # Test observation access
  first_obs <- mf$design$.obs[[1]]()
  expect_equal(dim(first_obs), c(3, 4))  # Check dimensions of first observation
  
  # Test with inconsistent dimensions
  x_bad <- list(
    matrix(1:12, 3, 4),
    matrix(13:24, 3, 4),
    matrix(25:36, 3, 4)   # Fixed dimensions to match
  )
  
  # This should now work since dimensions are consistent
  mf_consistent <- multiframe(x_bad, y)
  expect_s3_class(mf_consistent, "multiframe")
  
  # Test with different dimensions to ensure error
  x_bad2 <- list(
    matrix(1:12, 3, 4),
    matrix(13:24, 3, 4),
    matrix(25:40, 4, 4)   # Different number of rows
  )
  expect_error(
    multiframe(x_bad2, y),
    "All elements must have the same number of rows"
  )
  
  # Test with empty elements
  x_invalid <- list(matrix(1:12, 3, 4), numeric(0))
  expect_error(
    multiframe(x_invalid, y[1:2,]),
    "All elements must be non-empty"
  )
  
  # Test dimension mismatch
  expect_error(
    multiframe(x, y[1:2,]),
    "must be equal to"
  )
})

test_that("multiframe matrix constructor works correctly", {
  # Create test data
  x <- matrix(1:40, 10, 4)
  y <- data.frame(
    condition = rep(c("A", "B"), each=5),
    subject = rep(1:5, times=2)
  )
  
  # Test basic creation
  mf <- multiframe(x, y)
  expect_s3_class(mf, "multiframe")
  expect_equal(nrow(mf$design), 10)
  expect_true(all(c("condition", "subject", ".obs", ".index") %in% names(mf$design)))
  
  # Test observation access
  first_obs <- mf$design$.obs[[1]]()
  expect_equal(dim(first_obs), c(1, 4))  # Each observation is a 1x4 matrix
  
  # Test error conditions
  expect_error(
    multiframe(x[1:5,], y),
    "must be equal to"
  )
})

test_that("multiframe splitting works correctly", {
  # Create test data
  x <- matrix(1:40, 10, 4)
  y <- data.frame(
    condition = rep(c("A", "B"), each=5),
    subject = rep(1:5, times=2)
  )
  mf <- multiframe(x, y)
  
  # Test splitting by one variable
  split_cond <- split(mf, condition)
  expect_equal(nrow(split_cond), 2)  # Two conditions
  expect_equal(nrow(split_cond$data[[1]]), 5)  # 5 observations per condition
  
  # Test splitting by multiple variables
  split_both <- split(mf, condition, subject)
  expect_equal(nrow(split_both), 10)  # 2 conditions * 5 subjects
  
  # Test collapsing
  split_collapsed <- split(mf, condition, collapse=TRUE)
  expect_equal(nrow(split_collapsed), 2)
  expect_true(is.matrix(split_collapsed$data[[1]]))
})

test_that("multiframe summarization works correctly", {
  # Create test data
  x <- matrix(1:40, 10, 4)
  y <- data.frame(
    condition = rep(c("A", "B"), each=5),
    subject = rep(1:5, times=2)
  )
  mf <- multiframe(x, y)
  
  # Test basic summarization
  sum_cond <- summarize_by(mf, condition)
  expect_equal(nrow(sum_cond), 2)  # Two conditions
  expect_true(is.numeric(sum_cond$data[[1]]))
  expect_equal(length(sum_cond$data[[1]]), ncol(x))
  
  # Test with custom summary function
  sum_custom <- summarize_by(mf, condition, sfun=function(x) apply(x, 2, sd))
  expect_equal(nrow(sum_custom), 2)
  expect_equal(length(sum_custom$data[[1]]), ncol(x))
  
  # Test with data extraction
  sum_extract <- summarize_by(mf, condition, extract_data=TRUE)
  expect_true(is.matrix(sum_extract))
  expect_equal(nrow(sum_extract), 2)
  expect_equal(ncol(sum_extract), ncol(x))
})