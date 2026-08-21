make_boundary_hyperdesign <- function(first_cells = NULL, second_cells = NULL,
                                      space = "common") {
  column_design <- data.frame(axis = paste0("v", seq_len(4)))
  first <- multidesign(
    matrix(seq_len(12), nrow = 3),
    data.frame(entity = c("A", "B", "C")),
    column_design,
    cells = first_cells
  )
  second <- multidesign(
    matrix(seq_len(12) + 20L, nrow = 3),
    data.frame(entity = c("B", "C", "D")),
    column_design,
    cells = second_cells
  )
  hyperdesign(
    list(first = first, second = second),
    id = "entity",
    space = space
  )
}

test_that("as_multidesign row-stacks explicit and promoted masks", {
  mask <- matrix(seq_len(12) %% 4L != 0L, nrow = 3)
  hd <- make_boundary_hyperdesign(first_cells = mask)
  hd_before <- hd

  stacked <- as_multidesign(hd, .id = "source")
  expect_identical(
    cell_mask(stacked),
    rbind(mask, matrix(TRUE, nrow = 3, ncol = 4))
  )
  expect_equal(stacked$design$source, c(rep("first", 3), rep("second", 3)))
  expect_identical(hd, hd_before)

  unmasked <- make_boundary_hyperdesign()
  expect_null(cell_mask(as_multidesign(unmasked)))

  block_space <- make_boundary_hyperdesign(space = "block")
  expect_error(as_multidesign(block_space), "space = \\\"block\\\"")
})

test_that("init_transform fails closed or propagates output-shaped all-TRUE masks", {
  partial <- matrix(TRUE, nrow = 3, ncol = 4)
  partial[1, 1] <- FALSE
  hd_partial <- make_boundary_hyperdesign(first_cells = partial)
  expect_error(
    init_transform(hd_partial, multivarious::center()),
    "cannot transform partial cell masks"
  )

  all_true <- matrix(TRUE, nrow = 3, ncol = 4)
  hd_true <- make_boundary_hyperdesign(
    first_cells = all_true,
    second_cells = all_true
  )
  centered <- init_transform(hd_true, multivarious::center())
  expect_true(all(vapply(centered, has_cell_mask, logical(1))))
  expect_true(all(vapply(centered, function(block) all(cell_mask(block)), logical(1))))
  expect_true(all(vapply(
    centered,
    function(block) identical(dim(cell_mask(block)), dim(block$x)),
    logical(1)
  )))

  unmasked <- init_transform(make_boundary_hyperdesign(), multivarious::center())
  expect_true(all(vapply(unmasked, function(block) is.null(cell_mask(block)), logical(1))))
})
