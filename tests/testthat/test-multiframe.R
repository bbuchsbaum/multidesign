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

# --- Regression tests ---

test_that("obs_group matrix with custom ind validates using nrow", {
  X <- matrix(1:20, 5, 4)
  ind <- c(10, 20, 30, 40, 50)
  obs <- obs_group(X, ind=ind)
  expect_s3_class(obs, "observation_set")
  expect_length(obs, 5)
})

test_that("observation dispatches for numeric vectors", {
  v <- c(1.0, 2.0, 3.0)
  obs <- observation(v, 1)
  expect_true(is.function(obs))
  expect_equal(obs(), v)
})

test_that("observation.matrix returns proper matrix with drop=FALSE", {
  X <- matrix(1:20, 5, 4)
  obs <- observation(X, 3)
  result <- obs()
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 4)
  expect_equal(as.numeric(result), as.numeric(X[3, ]))
})

# --- Multiframe method parity tests ---

test_that("design.multiframe strips .obs and .index", {
  X <- matrix(rnorm(20 * 4), 20, 4)
  Y <- data.frame(condition = rep(c("A", "B"), each = 10), subject = rep(1:10, 2))
  mf <- multiframe(X, Y)
  des <- design(mf)
  expect_false(".obs" %in% names(des))
  expect_false(".index" %in% names(des))
  expect_true("condition" %in% names(des))
  expect_true("subject" %in% names(des))
  expect_equal(nrow(des), 20)
})

test_that("xdata.multiframe materializes correct matrix", {
  X <- matrix(1:40, 10, 4)
  Y <- data.frame(condition = rep(c("A", "B"), each = 5))
  mf <- multiframe(X, Y)
  mat <- xdata(mf)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(10, 4))
  expect_equal(mat, X)
})

test_that("split_indices.multiframe works", {
  X <- matrix(rnorm(40), 10, 4)
  Y <- data.frame(condition = rep(c("A", "B"), each = 5))
  mf <- multiframe(X, Y)
  si <- split_indices(mf, condition)
  expect_equal(nrow(si), 2)
  expect_true("indices" %in% names(si))
  expect_equal(sort(unlist(si$indices)), 1:10)
})

test_that("subset.multiframe filters correctly and returns NULL on empty", {
  X <- matrix(rnorm(40), 10, 4)
  Y <- data.frame(
    condition = rep(c("A", "B"), each = 5),
    subject = 1:10
  )
  mf <- multiframe(X, Y)

  mf_A <- subset(mf, condition == "A")
  expect_s3_class(mf_A, "multiframe")
  expect_equal(nrow(mf_A$design), 5)
  expect_true(all(mf_A$design$condition == "A"))

  # .index should be reset

  expect_equal(mf_A$design$.index, 1:5)

  # Observations should still work
  obs1 <- mf_A$design$.obs[[1]]()
  expect_true(is.matrix(obs1) || is.numeric(obs1))

  # Empty result
  mf_none <- subset(mf, condition == "C")
  expect_null(mf_none)
})

test_that("fold_over.multiframe creates valid folds with correct sizes", {
  X <- matrix(rnorm(40), 10, 4)
  Y <- data.frame(condition = rep(c("A", "B"), each = 5))
  mf <- multiframe(X, Y)

  folds <- fold_over(mf, condition)
  expect_s3_class(folds, "foldlist")
  expect_length(folds, 2)

  f1 <- folds[[1]]
  expect_s3_class(f1$analysis, "multiframe")
  expect_s3_class(f1$assessment, "multiframe")
  expect_equal(
    nrow(f1$analysis$design) + nrow(f1$assessment$design),
    10
  )
  # Each assessment set should have 5 observations
  expect_equal(nrow(f1$assessment$design), 5)
})

test_that("cv_rows.multiframe creates explicit row folds", {
  X <- matrix(1:24, 6, 4)
  Y <- data.frame(condition = letters[1:6])
  mf <- multiframe(X, Y)

  folds <- cv_rows(mf, rows = list(c(1, 2), c(5, 6)))
  expect_s3_class(folds, "foldlist")
  expect_length(folds, 2)

  f1 <- folds[[1]]
  expect_s3_class(f1$analysis, "multiframe")
  expect_s3_class(f1$assessment, "multiframe")
  expect_equal(xdata(f1$assessment), X[c(1, 2), , drop = FALSE])
  expect_equal(xdata(f1$analysis), X[-c(1, 2), , drop = FALSE])
  expect_equal(f1$held_out$rows, c(1L, 2L))
})

test_that("fold_over.multiframe can preserve original row ids", {
  X <- matrix(1:24, 6, 4)
  Y <- data.frame(condition = rep(c("A", "B"), each = 3))
  mf <- multiframe(X, Y)

  folds <- fold_over(mf, condition, preserve_row_ids = TRUE)
  f1 <- folds[[1]]

  expect_equal(f1$assessment$design$.orig_index, 1:3)
  expect_equal(f1$analysis$design$.orig_index, 4:6)
  expect_equal(f1$held_out$row_ids, 1:3)
})

test_that("[.observation_set extracts multiple observations", {
  X <- matrix(1:20, 5, 4)
  obs <- obs_group(X)

  # Extract subset using [ operator
  sub <- obs[1:3]
  expect_type(sub, "list")
  expect_length(sub, 3)
  # Each element should be an evaluated observation (matrix)
  expect_true(is.matrix(sub[[1]]))
  expect_equal(as.numeric(sub[[1]]), as.numeric(X[1, ]))
  expect_equal(as.numeric(sub[[3]]), as.numeric(X[3, ]))
})

test_that("observation.deflist creates lazy observations", {
  dl <- deflist::deflist(function(i) matrix(i * 10, 1, 3), len = 4)
  obs <- observation(dl, 2)
  expect_true(is.function(obs))
  result <- obs()
  expect_equal(result, dl[[2]])
})

test_that("print.observation_set does not error", {
  X <- matrix(rnorm(20), 5, 4)
  obs <- obs_group(X)
  expect_output(print(obs), "Observation Set")
})

test_that("print.multiframe does not error", {
  X <- matrix(rnorm(40), 10, 4)
  Y <- data.frame(
    condition = rep(c("A", "B"), each = 5),
    subject = 1:10
  )
  mf <- multiframe(X, Y)
  expect_output(print(mf), "Multiframe Object")
})
