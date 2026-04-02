library(tibble)

test_that("hyperdesign constructor works correctly", {
  # Create test multidesign objects
  X1 <- matrix(rnorm(50), 10, 5)
  Y1 <- tibble::tibble(
    condition = rep(c("A", "B"), each=5),
    subject = rep(1, 10),
    run = 1:10
  )
  col_design1 <- tibble::tibble(
    type = letters[1:5],
    group = rep(c("g1", "g2"), length.out=5)
  )
  d1 <- multidesign(X1, Y1, col_design1)

  X2 <- matrix(rnorm(50), 10, 5)
  Y2 <- tibble::tibble(
    condition = rep(c("A", "B"), each=5),
    subject = rep(2, 10),
    run = 1:10
  )
  d2 <- multidesign(X2, Y2, col_design1)  # Use the same column design

  # Create hyperdesign
  hd <- hyperdesign(list(d1, d2))

  # Test basic structure
  expect_s3_class(hd, "hyperdesign")
  expect_equal(length(hd), 2)
  expect_equal(attr(hd, "common_vars"), c("condition", "subject", "run"))

  # Test design extraction
  all_designs <- design(hd)
  expect_length(all_designs, 2)
  expect_equal(sort(names(all_designs[[1]])), sort(c(names(Y1))))

  # Test column_design extraction
  all_col_designs <- column_design(hd)
  expect_length(all_col_designs, 2)
  expect_equal(names(all_col_designs[[1]]), c("type", "group"))

  block1_col_design <- column_design(hd, block=1)
  expect_equal(names(block1_col_design), c("type", "group"))
  expect_equal(block1_col_design$type, letters[1:5])
})

test_that("fold_over creates valid cross-validation folds", {
  # Create test data with non-confounded variables
  X1 <- matrix(rnorm(50), 10, 5)
  Y1 <- tibble::tibble(
    condition = rep(c("A", "B"), each=5),
    subject = rep(c(1, 2), each=5),  # Not confounded with blocks
    run = 1:10
  )
  col_design1 <- tibble::tibble(
    type = letters[1:5],
    group = rep(c("g1", "g2"), length.out=5)
  )
  d1 <- multidesign(X1, Y1, col_design1)

  X2 <- matrix(rnorm(50), 10, 5)
  Y2 <- tibble::tibble(
    condition = rep(c("A", "B"), each=5),
    subject = rep(c(3, 4), each=5),  # Not confounded with blocks
    run = 1:10
  )
  d2 <- multidesign(X2, Y2, col_design1)

  # Create hyperdesign
  hd <- hyperdesign(list(d1, d2))

  # Test fold creation with non-confounded variable
  folds <- fold_over(hd, subject)
  expect_s3_class(folds, "foldlist")
  expect_length(folds, 4)  # 4 subjects = 4 folds

  # Test that each fold has both analysis and assessment sets
  expect_true(all(sapply(folds, function(f) {
    !is.null(f$analysis) && !is.null(f$assessment)
  })))

  # Test that each fold's assessment set contains only one subject
  expect_true(all(sapply(folds, function(f) {
    # The assessment set is a single multidesign object
    unique_subjects <- unique(f$assessment$design$subject)
    length(unique_subjects) == 1
  })))
})

test_that("fold_over handles confounded variables correctly", {
  # Create a simple design where 'group' is confounded with blocks
  X1 <- matrix(rnorm(20), 4, 5)
  Y1 <- tibble(
    group = c("A", "A", "A", "A"),  # Block 1 only has group A
    value = 1:4
  )

  X2 <- matrix(rnorm(20), 4, 5)
  Y2 <- tibble(
    group = c("B", "B", "B", "B"),  # Block 2 only has group B
    value = 5:8
  )

  md1 <- multidesign(X1, Y1)
  md2 <- multidesign(X2, Y2)
  hd <- hyperdesign(list(md1, md2))

  # Test error for confounded variable 'group'
  expect_error(
    fold_over(hd, group),
    "Variable 'group' is confounded with blocks"
  )

  # Test that non-confounded variable 'value' works
  expect_no_error(fold_over(hd, value))
})

test_that("cv_rows.hyperdesign supports synchronized row holdouts", {
  X1 <- matrix(1:30, 6, 5)
  Y1 <- tibble(run = 1:6, subject = "s1")
  X2 <- matrix(31:60, 6, 5)
  Y2 <- tibble(run = 1:6, subject = "s2")

  d1 <- multidesign(X1, Y1)
  d2 <- multidesign(X2, Y2)
  hd <- hyperdesign(list(d1, d2), block_names = c("subj1", "subj2"))

  folds <- cv_rows(hd, rows = list(
    list(subj1 = c(1, 2), subj2 = c(1, 2)),
    list(subj1 = c(3, 4), subj2 = c(3, 4))
  ))

  expect_s3_class(folds, "foldlist")
  expect_length(folds, 2)

  f1 <- folds[[1]]
  expect_s3_class(f1$analysis, "hyperdesign")
  expect_s3_class(f1$assessment, "hyperdesign")
  expect_equal(length(f1$assessment), 2)
  expect_equal(f1$assessment[[1]]$x, X1[c(1, 2), , drop = FALSE])
  expect_equal(f1$assessment[[2]]$x, X2[c(1, 2), , drop = FALSE])
  expect_equal(f1$analysis[[1]]$x, X1[-c(1, 2), , drop = FALSE])
  expect_equal(f1$analysis[[2]]$x, X2[-c(1, 2), , drop = FALSE])
  expect_equal(f1$held_out$subj1, c(1L, 2L))
  expect_equal(f1$held_out$subj2, c(1L, 2L))
})

test_that("preserve_row_ids carries source ids through hyperdesign folds", {
  X1 <- matrix(1:30, 6, 5)
  Y1 <- tibble(condition = rep(c("A", "B"), each = 3), subject = "s1")
  X2 <- matrix(31:60, 6, 5)
  Y2 <- tibble(condition = rep(c("A", "B"), each = 3), subject = "s2")

  d1 <- multidesign(X1, Y1)
  d2 <- multidesign(X2, Y2)
  hd <- hyperdesign(list(d1, d2), block_names = c("subj1", "subj2"))

  folds <- fold_over(hd, condition, preserve_row_ids = TRUE)
  f1 <- folds[[1]]
  expect_equal(f1$assessment$design$.orig_index, 1:3)
  expect_equal(f1$analysis[[1]]$design$.orig_index, 4:6)
  expect_equal(f1$held_out$row_ids, 1:3)

  row_folds <- cv_rows(
    hd,
    rows = list(list(subj1 = c(2, 3), subj2 = c(4, 5))),
    preserve_row_ids = TRUE
  )
  rf1 <- row_folds[[1]]
  expect_equal(rf1$assessment[[1]]$design$.orig_index, c(2L, 3L))
  expect_equal(rf1$assessment[[2]]$design$.orig_index, c(4L, 5L))
  expect_equal(rf1$held_out$row_ids$subj1, c(2L, 3L))
  expect_equal(rf1$held_out$row_ids$subj2, c(4L, 5L))
})

test_that("cv_rows.hyperdesign validates explicit row folds", {
  d1 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(run = 1:6))
  d2 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(run = 1:6))
  hd <- hyperdesign(list(d1, d2), block_names = c("subj1", "subj2"))

  expect_error(
    cv_rows(hd, rows = list(list(subj1 = c(1, 1)))),
    "must not contain duplicate row indices"
  )
  expect_error(
    cv_rows(hd, rows = list(list(subj1 = 1:6))),
    "cannot hold out all rows"
  )
  expect_error(
    cv_rows(hd, rows = list(list(unknown = 1:2))),
    "unknown hyperdesign blocks"
  )
})

test_that("data extraction methods work correctly", {
  # Create test data
  X1 <- matrix(1:20, 4, 5)
  Y1 <- tibble::tibble(
    condition = rep(c("A", "B"), each=2),
    subject = rep(1, 4),
    run = 1:4
  )
  col_design1 <- tibble::tibble(
    type = letters[1:5],
    group = rep(c("g1", "g2"), length.out=5)
  )
  d1 <- multidesign(X1, Y1, col_design1)

  X2 <- matrix(21:40, 4, 5)
  Y2 <- tibble::tibble(
    condition = rep(c("A", "B"), each=2),
    subject = rep(2, 4),
    run = 1:4
  )
  d2 <- multidesign(X2, Y2, col_design1)  # Use the same column design

  # Create hyperdesign
  hd <- hyperdesign(list(d1, d2))

  # Test xdata extraction
  all_data <- xdata(hd)
  expect_length(all_data, 2)
  expect_equal(all_data[[1]], X1)
  expect_equal(all_data[[2]], X2)

  block1_data <- xdata(hd, block=1)
  expect_equal(block1_data, X1)

  # Test design extraction
  all_designs <- design(hd)
  expect_length(all_designs, 2)
  expect_equal(sort(names(all_designs[[1]])), sort(c(names(Y1))))

  block1_design <- design(hd, block=1)
  expect_equal(sort(names(block1_design)), sort(c(names(Y1))))

  # Test column_design extraction
  all_col_designs <- column_design(hd)
  expect_length(all_col_designs, 2)
  expect_equal(names(all_col_designs[[1]]), c("type", "group"))

  block1_col_design <- column_design(hd, block=1)
  expect_equal(names(block1_col_design), c("type", "group"))
  expect_equal(block1_col_design$type, letters[1:5])
})

test_that("subsetting works correctly", {
  # Create test data
  X1 <- matrix(1:20, 4, 5)
  Y1 <- tibble::tibble(
    condition = rep(c("A", "B"), each=2),
    subject = rep(1, 4),
    block = c(1,1,2,2)
  )
  d1 <- multidesign(X1, Y1)

  X2 <- matrix(21:40, 4, 5)
  Y2 <- tibble::tibble(
    condition = rep(c("A", "B"), each=2),
    subject = rep(2, 4),
    block = c(1,1,2,2)
  )
  d2 <- multidesign(X2, Y2)

  hd <- hyperdesign(list(d1, d2))

  # Test subsetting by condition
  subset_A <- subset(hd, condition == "A")
  expect_s3_class(subset_A, "hyperdesign")
  expect_true(all(sapply(subset_A, function(x) all(x$design$condition == "A"))))

  # Test subsetting by multiple conditions
  subset_A1 <- subset(hd, condition == "A" & subject == 1)
  expect_true(all(sapply(subset_A1, function(x) {
    all(x$design$condition == "A" & x$design$subject == 1)
  })))
})

# --- Regression tests ---

test_that("print.hyperdesign does not error", {
  d1 <- multidesign(
    matrix(rnorm(50), 10, 5),
    data.frame(condition = rep(c("A", "B"), each=5))
  )
  d2 <- multidesign(
    matrix(rnorm(50), 10, 5),
    data.frame(condition = rep(c("A", "B"), each=5))
  )
  hd <- hyperdesign(list(d1, d2))
  expect_output(print(hd), "Hyperdesign Object")
})

test_that("df_to_hyperdesign creates hyperdesign from data frame", {
  sample_tibble <- tibble::tibble(
    felab = rep(1:2, each = 3),
    attention = rep(c("DA", "FA", "DA"), times = 2),
    subject = rep(1001:1002, each = 3),
    `v1` = rnorm(6),
    `v2` = rnorm(6),
    `v3` = rnorm(6)
  )

  hd <- df_to_hyperdesign(
    data = sample_tibble,
    design_vars = c("felab", "attention"),
    x_vars = c("v1", "v2", "v3"),
    split_var = "subject"
  )

  expect_s3_class(hd, "hyperdesign")
  expect_equal(length(hd), 2)
  # Each block should have 3 observations and 3 variables

  expect_equal(nrow(hd[[1]]$x), 3)
  expect_equal(ncol(hd[[1]]$x), 3)
  # Design should have felab and attention columns
  expect_true(all(c("felab", "attention") %in% names(hd[[1]]$design)))
})

test_that("block_indices.hyperdesign works with and without i", {
  d1 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(cond = rep(c("A","B"), 3)))
  d2 <- multidesign(matrix(rnorm(40), 8, 5), data.frame(cond = rep(c("A","B"), 4)))
  hd <- hyperdesign(list(d1, d2))

  # Without i: returns list of all block indices
  all_idx <- block_indices(hd)
  expect_type(all_idx, "list")
  expect_length(all_idx, 2)
  expect_equal(all_idx[[1]], 1:5)
  expect_equal(all_idx[[2]], 6:10)

  # With i: returns indices for specific block
  idx1 <- block_indices(hd, 1)
  expect_equal(idx1, 1:5)

  # byrow = TRUE
  all_row_idx <- block_indices(hd, byrow = TRUE)
  expect_type(all_row_idx, "list")
  expect_length(all_row_idx, 2)
  expect_equal(all_row_idx[[1]], 1:6)
  expect_equal(all_row_idx[[2]], 7:14)

  row_idx1 <- block_indices(hd, 1, byrow = TRUE)
  expect_equal(row_idx1, 1:6)
})

test_that("fold_over.hyperdesign leave-one-block-out works", {
  d1 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(cond = rep(c("A","B"), 3)))
  d2 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(cond = rep(c("A","B"), 3)))
  d3 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(cond = rep(c("A","B"), 3)))
  hd <- hyperdesign(list(d1, d2, d3))

  # Leave-one-block-out (no split vars)
  folds <- fold_over(hd)
  expect_s3_class(folds, "foldlist")
  expect_length(folds, 3)

  f1 <- folds[[1]]
  # Assessment is a multidesign (single block)
  expect_s3_class(f1$assessment, "multidesign")
  # Analysis is a hyperdesign (remaining blocks)
  expect_s3_class(f1$analysis, "hyperdesign")
  expect_equal(length(f1$analysis), 2)
  expect_equal(nrow(f1$assessment$x), 6)
})

test_that("print.foldlist does not error", {
  d1 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(cond = rep(c("A","B"), 3)))
  d2 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(cond = rep(c("A","B"), 3)))
  hd <- hyperdesign(list(d1, d2))

  folds <- fold_over(hd)
  expect_output(print(folds), "Cross-Validation Folds")

  # Also test with variable-based folds
  folds2 <- fold_over(hd, cond)
  expect_output(print(folds2), "Cross-Validation Folds")
})

test_that("print.foldlist formats nested held-out metadata", {
  d1 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(run = 1:6))
  d2 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(run = 1:6))
  hd <- hyperdesign(list(d1, d2), block_names = c("subj1", "subj2"))

  folds <- cv_rows(hd, rows = list(list(subj1 = 1:2, subj2 = 1:2)))
  expect_output(print(folds), "subj1")
  expect_output(print(folds), "subj2")
})

test_that("subset.hyperdesign errors when nothing matches", {
  d1 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(cond = rep(c("A","B"), 3)))
  d2 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(cond = rep(c("A","B"), 3)))
  hd <- hyperdesign(list(d1, d2))

  expect_error(subset(hd, cond == "Z"), "does not match any rows")
})
