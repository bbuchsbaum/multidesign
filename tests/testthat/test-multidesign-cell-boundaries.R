make_boundary_multidesign <- function(offset = 0L, cells = NULL) {
  x <- matrix(seq_len(12) + offset, nrow = 3, ncol = 4)
  multidesign(
    x,
    data.frame(group = rep(if (offset == 0L) "A" else "B", 3)),
    data.frame(variable = paste0("v", seq_len(4))),
    cells = cells
  )
}

test_that("bind_multidesign promotes only mixed NULL masks", {
  mask <- matrix(seq_len(12) %% 3L != 0L, nrow = 3)
  masked <- make_boundary_multidesign(cells = mask)
  unmasked <- make_boundary_multidesign(20L)
  masked_before <- masked
  unmasked_before <- unmasked

  mixed <- bind_multidesign(masked, unmasked, .id = "source")
  expect_identical(
    cell_mask(mixed),
    rbind(mask, matrix(TRUE, nrow = 3, ncol = 4))
  )
  expect_equal(mixed$x, rbind(masked$x, unmasked$x))
  expect_equal(mixed$design$source, c(rep(1L, 3), rep(2L, 3)))

  all_null <- bind_multidesign(unmasked, unmasked)
  expect_null(cell_mask(all_null))
  expect_null(cell_mask(bind_multidesign(list(unmasked))))

  expect_identical(masked, masked_before)
  expect_identical(unmasked, unmasked_before)
})

test_that("reduce fails closed for partial masks and propagates all-TRUE masks", {
  partial <- make_boundary_multidesign(
    cells = matrix(c(rep(TRUE, 11), FALSE), nrow = 3)
  )
  expect_error(reduce.multidesign(partial, nc = 2), "cannot transform partial cell masks")

  all_true <- make_boundary_multidesign(cells = matrix(TRUE, nrow = 3, ncol = 4))
  reduced <- reduce.multidesign(all_true, nc = 2)
  expect_s3_class(reduced, "reduced_multidesign")
  expect_identical(dim(cell_mask(reduced)), dim(reduced$x))
  expect_true(all(cell_mask(reduced)))
  expect_equal(nrow(reduced$column_design), ncol(reduced$x))

  unmasked <- reduce.multidesign(make_boundary_multidesign(), nc = 2)
  expect_null(cell_mask(unmasked))
  expect_equal(nrow(unmasked$column_design), ncol(unmasked$x))
})
