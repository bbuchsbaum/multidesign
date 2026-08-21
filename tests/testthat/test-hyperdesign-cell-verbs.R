make_hyper_cell_block <- function(ids, keep, offset = 0L) {
  x <- matrix(seq_len(length(ids) * 3L) + offset, nrow = length(ids), ncol = 3L)
  cells <- matrix((seq_len(length(x)) + offset) %% 3L != 0L, nrow = length(ids))
  multidesign(
    x,
    data.frame(entity = ids, keep = keep),
    data.frame(axis = c("x", "y", "z")),
    cells = cells
  )
}

test_that("hyperdesign subset preserves block masks and contracts", {
  one <- make_hyper_cell_block(c("A", "B", "C"), c(TRUE, FALSE, TRUE))
  two <- make_hyper_cell_block(c("B", "C", "D"), c(TRUE, TRUE, FALSE), 20L)
  hd <- hyperdesign(list(one = one, two = two), id = "entity", space = "common")

  filtered <- subset(hd, keep)

  expect_identical(cell_mask(filtered$one), cell_mask(one)[c(1, 3), , drop = FALSE])
  expect_identical(cell_mask(filtered$two), cell_mask(two)[c(1, 2), , drop = FALSE])
  expect_equal(entity_id(filtered), "entity")
  expect_equal(column_space(filtered), "common")
  expect_equal(correspondence(filtered)$global_ids, c("A", "C", "B"))
})

test_that("hyperdesign entity folds and cv_rows preserve per-block masks", {
  one <- make_hyper_cell_block(c("A", "B", "C"), rep(TRUE, 3))
  two <- make_hyper_cell_block(c("B", "C", "D"), rep(TRUE, 3), 20L)
  hd <- hyperdesign(list(one = one, two = two), id = "entity", space = "common")

  entity_fold <- fold_over(hd, entity)[[2]]
  expect_equal(entity_fold$held_out$entity, "B")
  expect_identical(cell_mask(entity_fold$assessment$one), cell_mask(one)[2, , drop = FALSE])
  expect_identical(cell_mask(entity_fold$assessment$two), cell_mask(two)[1, , drop = FALSE])
  expect_identical(cell_mask(entity_fold$analysis$one), cell_mask(one)[c(1, 3), , drop = FALSE])
  expect_identical(cell_mask(entity_fold$analysis$two), cell_mask(two)[c(2, 3), , drop = FALSE])

  explicit <- cv_rows(hd, list(list(one = c(1, 3), two = 2L)))[[1]]
  expect_identical(cell_mask(explicit$assessment$one), cell_mask(one)[c(1, 3), , drop = FALSE])
  expect_identical(cell_mask(explicit$assessment$two), cell_mask(two)[2, , drop = FALSE])
  expect_identical(cell_mask(explicit$analysis$one), cell_mask(one)[2, , drop = FALSE])
  expect_identical(cell_mask(explicit$analysis$two), cell_mask(two)[c(1, 3), , drop = FALSE])
})
