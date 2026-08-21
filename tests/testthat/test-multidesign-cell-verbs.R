make_cell_verb_design <- function() {
  x <- matrix(seq_len(30), nrow = 6, ncol = 5)
  design <- data.frame(
    group = rep(c("A", "B"), each = 3),
    keep = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  )
  column_design <- data.frame(
    variable = paste0("v", seq_len(5)),
    region = c("left", "right", "left", "right", "left")
  )
  cells <- matrix(
    (seq_len(length(x)) %% 4L) != 0L,
    nrow = nrow(x),
    ncol = ncol(x)
  )
  cells[3, ] <- FALSE
  multidesign(x, design, column_design, cells = cells)
}

test_that("row and column verbs slice cells with the data", {
  md <- make_cell_verb_design()
  mask <- cell_mask(md)

  filtered <- subset(md, keep)
  expect_identical(cell_mask(filtered), mask[c(1, 3, 5, 6), , drop = FALSE])
  expect_false(any(cell_mask(filtered)[2, ]))

  groups <- split(md, group)
  expect_identical(cell_mask(groups[[1]]), mask[1:3, , drop = FALSE])
  expect_identical(cell_mask(groups[[2]]), mask[4:6, , drop = FALSE])

  selected <- select_variables(md, region == "left")
  expect_identical(cell_mask(selected), mask[, c(1, 3, 5), drop = FALSE])

  composed <- select_variables(filtered, region == "left")
  expect_identical(
    cell_mask(composed),
    mask[c(1, 3, 5, 6), c(1, 3, 5), drop = FALSE]
  )
})

test_that("fold and explicit-row CV preserve exact analysis and assessment masks", {
  md <- make_cell_verb_design()
  mask <- cell_mask(md)

  row_folds <- cv_rows(
    md,
    rows = list(c(2, 5)),
    preserve_row_ids = TRUE
  )
  first <- row_folds[[1]]
  second_materialization <- row_folds[[1]]
  expect_identical(cell_mask(first$assessment), mask[c(2, 5), , drop = FALSE])
  expect_identical(cell_mask(first$analysis), mask[c(1, 3, 4, 6), , drop = FALSE])
  expect_identical(cell_mask(second_materialization$assessment), cell_mask(first$assessment))
  expect_equal(first$assessment$design$.orig_index, c(2L, 5L))

  group_folds <- fold_over(md, group)
  expect_identical(cell_mask(group_folds[[1]]$assessment), mask[1:3, , drop = FALSE])
  expect_identical(cell_mask(group_folds[[1]]$analysis), mask[4:6, , drop = FALSE])
})

test_that("NULL masks remain NULL through supported verbs", {
  source <- make_cell_verb_design()
  md <- multidesign(source$x, design(source), source$column_design)

  expect_null(cell_mask(subset(md, keep)))
  expect_true(all(vapply(split(md, group), function(part) is.null(cell_mask(part)), logical(1))))
  expect_null(cell_mask(select_variables(md, region == "left")))
  expect_null(cell_mask(cv_rows(md, list(1:2))[[1]]$assessment))
  expect_null(cell_mask(fold_over(md, group)[[1]]$analysis))
})

test_that("random valid selections agree with direct mask indexing", {
  md <- make_cell_verb_design()
  mask <- cell_mask(md)
  set.seed(20260821)

  for (iteration in seq_len(30L)) {
    rows <- sort(sample(seq_len(nrow(md$x)), sample(seq_len(nrow(md$x)), 1L)))
    cols <- sort(sample(seq_len(ncol(md$x)), sample(seq_len(ncol(md$x)), 1L)))

    row_design <- design(md)
    row_design$chosen <- seq_len(nrow(row_design)) %in% rows
    candidate <- multidesign(md$x, row_design, md$column_design, cells = mask)
    row_selected <- subset(candidate, chosen)

    col_design <- row_selected$column_design
    col_design$chosen <- seq_len(nrow(col_design)) %in% cols
    row_selected$column_design <- tibble::as_tibble(col_design)
    result <- select_variables(row_selected, chosen)

    expect_identical(cell_mask(result), mask[rows, cols, drop = FALSE])
  }
})
