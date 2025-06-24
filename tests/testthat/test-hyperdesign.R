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
