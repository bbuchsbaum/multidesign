library(testthat)
library(tibble)

test_that("cross_validate basic pipeline with multidesign folds", {
  set.seed(42)
  X <- matrix(rnorm(100 * 5), 100, 5)
  Y <- data.frame(group = rep(c("A", "B"), each = 50))
  mds <- multidesign(X, Y)
  folds <- fold_over(mds, group)

  result <- cross_validate(
    folds,
    fit_fn = function(analysis) {
      list(mean = colMeans(xdata(analysis)))
    },
    score_fn = function(model, assessment) {
      pred_error <- mean((xdata(assessment) -
        matrix(model$mean,
          nrow = nrow(xdata(assessment)),
          ncol = length(model$mean), byrow = TRUE
        ))^2)
      c(mse = pred_error)
    }
  )

  expect_s3_class(result, "cv_result")
  expect_true(".fold" %in% names(result$scores))
  expect_true("mse" %in% names(result$scores))
  expect_equal(nrow(result$scores), 2)
  expect_true(all(result$scores$mse > 0))
})

test_that("cross_validate with named multi-score return", {
  set.seed(42)
  X <- matrix(rnorm(60 * 4), 60, 4)
  Y <- data.frame(cond = rep(c("A", "B", "C"), each = 20))
  mds <- multidesign(X, Y)
  folds <- fold_over(mds, cond)

  result <- cross_validate(
    folds,
    fit_fn = function(analysis) list(m = colMeans(xdata(analysis))),
    score_fn = function(model, assessment) {
      mat <- xdata(assessment)
      c(
        mse = mean((mat - matrix(model$m, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE))^2),
        mae = mean(abs(mat - matrix(model$m, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)))
      )
    }
  )

  expect_s3_class(result, "cv_result")
  expect_true("mse" %in% names(result$scores))
  expect_true("mae" %in% names(result$scores))
  expect_equal(nrow(result$scores), 3)
})

test_that("cross_validate with hyperdesign leave-one-block-out folds", {
  set.seed(42)
  d1 <- multidesign(
    matrix(rnorm(20), 10, 2),
    data.frame(condition = rep(c("A", "B"), 5))
  )
  d2 <- multidesign(
    matrix(rnorm(20), 10, 2),
    data.frame(condition = rep(c("A", "B"), 5))
  )
  hd <- hyperdesign(list(d1, d2))
  folds <- fold_over(hd)

  result <- cross_validate(
    folds,
    fit_fn = function(analysis) {
      # analysis is a hyperdesign; extract all data
      all_x <- do.call(rbind, lapply(analysis, function(b) xdata(b)))
      list(mean = colMeans(all_x))
    },
    score_fn = function(model, assessment) {
      mat <- xdata(assessment)
      c(mse = mean((mat - matrix(model$mean, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE))^2))
    }
  )

  expect_s3_class(result, "cv_result")
  expect_equal(nrow(result$scores), 2)
})

test_that("cross_validate errors on bad score_fn return type", {
  set.seed(42)
  X <- matrix(rnorm(40), 10, 4)
  Y <- data.frame(group = rep(c("A", "B"), each = 5))
  mds <- multidesign(X, Y)
  folds <- fold_over(mds, group)

  expect_error(
    cross_validate(
      folds,
      fit_fn = function(analysis) list(),
      score_fn = function(model, assessment) "not_numeric"
    ),
    "score_fn must return"
  )
})

test_that("print.cv_result works", {
  scores <- tibble(.fold = 1:3, mse = c(0.5, 0.6, 0.55))
  cr <- new_cv_result(scores)
  expect_output(print(cr), "Cross-Validation Result")
  expect_output(print(cr), "Folds:")
  expect_output(print(cr), "mse")
})

test_that("summary.cv_result returns correct tibble", {
  scores <- tibble(.fold = 1:4, mse = c(1.0, 2.0, 3.0, 4.0), acc = c(0.9, 0.8, 0.85, 0.95))
  cr <- new_cv_result(scores)
  s <- summary(cr)
  expect_s3_class(s, "tbl_df")
  expect_equal(nrow(s), 2)
  expect_true(all(c("metric", "mean", "sd", "min", "max", "median") %in% names(s)))

  mse_row <- s[s$metric == "mse", ]
  expect_equal(mse_row$mean, 2.5)
  expect_equal(mse_row$min, 1.0)
  expect_equal(mse_row$max, 4.0)
  expect_equal(mse_row$median, 2.5)
})

test_that("cross_validate with multiframe folds", {
  set.seed(42)
  X <- matrix(rnorm(40), 10, 4)
  Y <- data.frame(condition = rep(c("A", "B"), each = 5))
  mf <- multiframe(X, Y)
  folds <- fold_over(mf, condition)

  result <- cross_validate(
    folds,
    fit_fn = function(analysis) {
      list(mean = colMeans(xdata(analysis)))
    },
    score_fn = function(model, assessment) {
      mat <- xdata(assessment)
      c(mse = mean((mat - matrix(model$mean, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE))^2))
    }
  )

  expect_s3_class(result, "cv_result")
  expect_equal(nrow(result$scores), 2)
  expect_true("mse" %in% names(result$scores))
})
