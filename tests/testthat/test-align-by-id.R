make_alignment_block <- function(ids, values, cells = NULL) {
  axes <- colnames(values)
  if (is.null(axes)) {
    axes <- paste0("v", seq_len(ncol(values)))
  }
  multidesign(
    values,
    data.frame(entity = ids),
    data.frame(axis = axes),
    cells = cells
  )
}

reference_alignment <- function(hd, fill = NA_real_, sort_ids = FALSE) {
  cr <- correspondence(hd, sort_ids = sort_ids)
  output <- array(fill, c(cr$n_global, ncol(hd[[1]]$x), length(hd)))
  output_cells <- array(FALSE, dim(output))
  for (block_index in seq_along(hd)) {
    block <- hd[[block_index]]
    mask <- cell_mask(block)
    if (is.null(mask)) {
      mask <- matrix(TRUE, nrow(block$x), ncol(block$x))
    }
    for (local_row in seq_len(nrow(block$x))) {
      global_row <- cr$row_map[[block_index]][[local_row]]
      for (coordinate in seq_len(ncol(block$x))) {
        if (mask[local_row, coordinate]) {
          output[global_row, coordinate, block_index] <- block$x[local_row, coordinate]
          output_cells[global_row, coordinate, block_index] <- TRUE
        }
      }
    }
  }
  list(x = output, cells = output_cells, observed = cr$overlap$incidence, correspondence = cr)
}

test_that("align_by_id separates row incidence, cell masks, and observed NA", {
  first_x <- matrix(
    c(1, 2, NA, 4, 5, 6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("x", "y"))
  )
  first_cells <- matrix(
    c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
    nrow = 3,
    byrow = TRUE
  )
  second_x <- matrix(
    c(10, 11, 12, 13, 14, 15),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("x", "y"))
  )
  first <- make_alignment_block(c("B", "A", "C"), first_x, first_cells)
  second <- make_alignment_block(c("B", "C", "D"), second_x)
  hd <- hyperdesign(list(first = first, second = second), id = "entity", space = "common")
  before <- hd

  aligned <- align_by_id(hd)

  expect_s3_class(aligned, "aligned_hyperdesign")
  expect_identical(dim(aligned$x), c(4L, 2L, 2L))
  expect_equal(dimnames(aligned$x)$entity, c("B", "A", "C", "D"))
  expect_equal(dimnames(aligned$x)$variable, c("x", "y"))
  expect_equal(dimnames(aligned$x)$block, c("first", "second"))
  expect_identical(unname(aligned$observed[, "first"]), c(TRUE, TRUE, TRUE, FALSE))
  expect_identical(unname(aligned$observed[, "second"]), c(TRUE, FALSE, TRUE, TRUE))

  expect_equal(aligned$x["B", "x", "first"], 1)
  expect_true(is.na(aligned$x["B", "y", "first"]))
  expect_false(aligned$cells["B", "y", "first"])
  expect_true(is.na(aligned$x["A", "x", "first"]))
  expect_true(aligned$cells["A", "x", "first"])
  expect_true(aligned$observed["C", "first"])
  expect_false(any(aligned$cells["C", , "first"]))
  expect_false(aligned$observed["D", "first"])
  expect_true(all(is.na(aligned$x["D", , "first"])))
  expect_identical(hd, before)

  printed <- testthat::capture_output(print(aligned))
  expect_match(paste(printed, collapse = "\n"), "4 entities x 2 variables x 2 blocks")
})

test_that("custom fill is presentation only", {
  first_x <- matrix(c(NA, 2, 3, 4), nrow = 2, byrow = TRUE)
  colnames(first_x) <- c("x", "y")
  cells <- matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2, byrow = TRUE)
  first <- make_alignment_block(c("A", "B"), first_x, cells)
  second <- make_alignment_block(c("B"), matrix(c(5, 6), nrow = 1, dimnames = list(NULL, c("x", "y"))))
  hd <- hyperdesign(list(first = first, second = second), id = "entity", space = "common")

  aligned <- align_by_id(hd, fill = -999)
  expect_true(is.na(aligned$x["A", "x", "first"]))
  expect_true(aligned$cells["A", "x", "first"])
  expect_equal(aligned$x["A", "y", "first"], -999)
  expect_false(aligned$cells["A", "y", "first"])
  expect_equal(aligned$x["A", "x", "second"], -999)
  expect_false(aligned$observed["A", "second"])
  expect_error(align_by_id(hd, fill = c(0, 1)), "one numeric")
  expect_error(align_by_id(hd, fill = "missing"), "one numeric")
})

test_that("alignment validates correspondence and common space", {
  values <- matrix(seq_len(4), nrow = 2)
  first <- make_alignment_block(c("A", "B"), values)
  second <- make_alignment_block(c("A", "C"), values + 10)

  expect_error(
    align_by_id(hyperdesign(list(first = first, second = second), space = "common")),
    "requires declared correspondence"
  )
  expect_error(
    align_by_id(hyperdesign(list(first = first, second = second), id = "entity")),
    'space = "common"'
  )
  expect_error(
    align_by_id(hyperdesign(list(first = first, second = second), id = "entity", space = "block")),
    'space = "common"'
  )
})

test_that("sorted, positional, one-block, and disconnected alignment are supported", {
  column_design <- data.frame(axis = c("x", "y"))
  one <- multidesign(matrix(1:4, nrow = 2), data.frame(entity = c("D", "B")), column_design)
  two <- multidesign(matrix(5:8, nrow = 2), data.frame(entity = c("A", "C")), column_design)
  disconnected <- hyperdesign(list(one = one, two = two), id = "entity", space = "common")
  aligned <- align_by_id(disconnected, sort_ids = TRUE)
  expect_equal(dimnames(aligned$x)$entity, c("A", "B", "C", "D"))
  expect_false(aligned$correspondence$overlap$connected)

  positional <- hyperdesign(
    list(one = one, two = one),
    positional = TRUE,
    space = "common"
  )
  expect_equal(dim(align_by_id(positional)$x), c(2L, 2L, 2L))

  single <- hyperdesign(list(one = one), id = "entity", space = "common")
  single_aligned <- align_by_id(single)
  expect_equal(dim(single_aligned$x), c(2L, 2L, 1L))
  expect_true(all(single_aligned$observed))
})

test_that("aligned arrays reconstruct every local observed value", {
  set.seed(20260823)
  for (iteration in seq_len(25L)) {
    universe <- paste0("E", seq_len(sample(5:10, 1L)))
    n_blocks <- sample(2:5, 1L)
    column_design <- data.frame(axis = paste0("v", seq_len(3L)))
    blocks <- lapply(seq_len(n_blocks), function(block_index) {
      ids <- sample(universe, sample(seq_along(universe), 1L))
      values <- matrix(rnorm(length(ids) * 3L), nrow = length(ids), ncol = 3L)
      cells <- matrix(sample(c(TRUE, FALSE), length(values), replace = TRUE), nrow = length(ids))
      multidesign(values, data.frame(entity = ids), column_design, cells = cells)
    })
    names(blocks) <- paste0("block", seq_along(blocks))
    hd <- hyperdesign(blocks, id = "entity", space = "common")

    actual <- align_by_id(hd, fill = -7, sort_ids = TRUE)
    reference <- reference_alignment(hd, fill = -7, sort_ids = TRUE)
    expect_equal(unname(actual$x), reference$x, tolerance = 0)
    expect_identical(unname(actual$cells), reference$cells)
    expect_identical(actual$observed, reference$observed)

    for (block_index in seq_along(hd)) {
      rows <- actual$correspondence$row_map[[block_index]]
      local_mask <- cell_mask(hd[[block_index]])
      local_values <- actual$x[rows, , block_index, drop = FALSE]
      local_cells <- actual$cells[rows, , block_index, drop = FALSE]
      dim(local_values) <- dim(hd[[block_index]]$x)
      dim(local_cells) <- dim(local_mask)
      expect_identical(local_cells, unname(local_mask))
      expect_equal(
        unname(local_values[local_mask]),
        unname(hd[[block_index]]$x[local_mask]),
        tolerance = 0
      )
    }
  }
})
