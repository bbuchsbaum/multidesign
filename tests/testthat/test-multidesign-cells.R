test_that("cell masks are explicit and backward compatible", {
  x <- matrix(seq_len(12), nrow = 3)
  design <- data.frame(entity = letters[1:3])
  cells <- matrix(
    c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
    nrow = 3
  )

  legacy <- multidesign(x, design)
  masked <- multidesign(x, design, cells = cells)

  expect_false(has_cell_mask(legacy))
  expect_null(cell_mask(legacy))
  expect_false("cells" %in% names(legacy))
  expect_true(has_cell_mask(masked))
  expect_identical(cell_mask(masked), cells)

  old_shape <- legacy
  old_shape$cells <- NULL
  expect_null(cell_mask(old_shape))
  expect_false(has_cell_mask(old_shape))
})

test_that("cell-mask validation rejects ambiguous representations", {
  x <- matrix(seq_len(12), nrow = 3)
  design <- data.frame(entity = letters[1:3])

  expect_error(multidesign(x, design, cells = rep(TRUE, length(x))), "logical matrix")
  expect_error(multidesign(x, design, cells = matrix(1, 3, 4)), "logical matrix")
  expect_error(multidesign(x, design, cells = array(TRUE, c(3, 4, 1))), "logical matrix")
  expect_error(multidesign(x, design, cells = matrix(TRUE, 2, 6)), "dimensions identical")

  with_na <- matrix(TRUE, nrow = 3, ncol = 4)
  with_na[2, 3] <- NA
  expect_error(multidesign(x, design, cells = with_na), "cannot contain NA")
})

test_that("all-FALSE rows remain and x values never define the mask", {
  x <- matrix(seq_len(12), nrow = 3)
  x[1, 1] <- NA_real_
  x[2, 2] <- NaN
  x[3, 3] <- Inf
  cells <- matrix(TRUE, nrow = 3, ncol = 4)
  cells[2, ] <- FALSE
  design <- data.frame(entity = letters[1:3])

  masked <- multidesign(x, design, cells = cells)

  expect_equal(nrow(masked$x), 3L)
  expect_identical(cell_mask(masked), cells)
  expect_true(cell_mask(masked)[1, 1])
  expect_true(cell_mask(masked)[3, 3])
  expect_false(any(cell_mask(masked)[2, ]))
})

test_that("printing reports masks only when explicitly present", {
  x <- matrix(seq_len(6), nrow = 2)
  design <- data.frame(entity = c("A", "B"))
  cells <- matrix(c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE), nrow = 2)

  legacy_output <- crayon::strip_style(paste(
    testthat::capture_output(print(multidesign(x, design))),
    collapse = "\n"
  ))
  masked_output <- crayon::strip_style(paste(
    testthat::capture_output(print(multidesign(x, design, cells = cells))),
    collapse = "\n"
  ))

  expect_false(grepl("Observed Cells", legacy_output, fixed = TRUE))
  expect_match(masked_output, "Observed Cells")
  expect_match(masked_output, "3 of 6 explicitly observed")

  reduced <- multidesign(x, design, cells = cells)
  class(reduced) <- c("reduced_multidesign", "multidesign")
  reduced$projector <- structure(list(), class = "projector")
  expect_true(has_cell_mask(reduced))
})
