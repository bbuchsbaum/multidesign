library(testthat)
library(tibble)

test_that("select_variables.multidesign filters columns by single condition", {
  X <- matrix(rnorm(20 * 10), 20, 10)
  Y <- tibble(condition = rep(c("A", "B"), each = 10))
  col_info <- tibble(
    region = rep(c("frontal", "parietal"), 5),
    hemisphere = rep(c("left", "right"), each = 5)
  )
  mds <- multidesign(X, Y, col_info)

  mds_frontal <- select_variables(mds, region == "frontal")
  expect_s3_class(mds_frontal, "multidesign")
  expect_equal(ncol(mds_frontal$x), 5)
  expect_true(all(mds_frontal$column_design$region == "frontal"))
  expect_equal(nrow(mds_frontal$x), 20)
})

test_that("select_variables.multidesign filters columns by multiple conditions", {
  X <- matrix(rnorm(20 * 10), 20, 10)
  Y <- tibble(condition = rep(c("A", "B"), each = 10))
  col_info <- tibble(
    region = rep(c("frontal", "parietal"), 5),
    hemisphere = rep(c("left", "right"), each = 5)
  )
  mds <- multidesign(X, Y, col_info)

  mds_sub <- select_variables(mds, region == "frontal" & hemisphere == "left")
  expect_s3_class(mds_sub, "multidesign")
  expect_true(all(mds_sub$column_design$region == "frontal"))
  expect_true(all(mds_sub$column_design$hemisphere == "left"))
  expect_true(ncol(mds_sub$x) > 0)
})

test_that("select_variables.multidesign errors when no columns match", {
  X <- matrix(rnorm(20 * 10), 20, 10)
  Y <- tibble(condition = rep(c("A", "B"), each = 10))
  col_info <- tibble(
    region = rep(c("frontal", "parietal"), 5)
  )
  mds <- multidesign(X, Y, col_info)

  expect_error(
    select_variables(mds, region == "occipital"),
    "no columns match"
  )
})

test_that("select_variables.multidesign selects correct column data", {
  X <- matrix(1:30, 3, 10)
  Y <- tibble(row_id = 1:3)
  col_info <- tibble(
    var_name = paste0("v", 1:10),
    group = rep(c("A", "B"), each = 5)
  )
  mds <- multidesign(X, Y, col_info)

  mds_A <- select_variables(mds, group == "A")
  expect_equal(ncol(mds_A$x), 5)
  expect_equal(mds_A$x, X[, 1:5, drop = FALSE])
  expect_equal(mds_A$column_design$var_name, paste0("v", 1:5))
})

test_that("select_variables.hyperdesign applies per-block", {
  col_info <- tibble(
    region = rep(c("frontal", "parietal"), 5),
    hemisphere = rep(c("left", "right"), each = 5)
  )
  d1 <- multidesign(
    matrix(rnorm(10 * 10), 10, 10),
    data.frame(condition = rep(c("A", "B"), 5)),
    col_info
  )
  d2 <- multidesign(
    matrix(rnorm(10 * 10), 10, 10),
    data.frame(condition = rep(c("A", "B"), 5)),
    col_info
  )
  hd <- hyperdesign(list(d1, d2))

  hd_frontal <- select_variables(hd, region == "frontal")
  expect_s3_class(hd_frontal, "hyperdesign")
  expect_equal(length(hd_frontal), 2)
  expect_true(all(sapply(hd_frontal, function(b) ncol(b$x)) == 5))
  expect_true(all(sapply(hd_frontal, function(b) all(b$column_design$region == "frontal"))))
})
