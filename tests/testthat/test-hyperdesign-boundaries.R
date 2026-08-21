make_boundary_block <- function(ids, ncol = 2L, offset = 0L,
                                column_design = NULL) {
  if (is.null(column_design)) {
    column_design <- data.frame(axis = seq_len(ncol))
  }
  multidesign(
    matrix(seq_len(length(ids) * ncol) + offset, nrow = length(ids)),
    data.frame(entity = ids, condition = rep("A", length(ids))),
    column_design
  )
}

test_that("as_multidesign keeps stacking semantics for unset and common space", {
  d1 <- make_boundary_block(c("A", "B"))
  d2 <- make_boundary_block(c("A", "C"), offset = 20L)

  legacy <- as_multidesign(hyperdesign(list(one = d1, two = d2)))
  common <- as_multidesign(
    hyperdesign(list(one = d1, two = d2), id = "entity", space = "common"),
    .id = "source"
  )

  expect_equal(nrow(legacy$x), 4L)
  expect_equal(nrow(common$x), 4L)
  expect_equal(common$design$entity, c("A", "B", "A", "C"))
  expect_equal(common$design$source, c("one", "one", "two", "two"))
})

test_that("space block refuses every hyperdesign stacking path", {
  d1 <- make_boundary_block(c("A", "B"))
  d2 <- make_boundary_block(c("A", "C"), offset = 20L)
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity", space = "block")

  expect_error(as_multidesign(hd), 'space = "block"', fixed = TRUE)
  expect_error(bind_multidesign(hd), 'space = "block"', fixed = TRUE)
  expect_error(bind_multidesign(hd, d1), 'space = "block"', fixed = TRUE)
})

test_that("bind_multidesign still accepts one multidesign or a plain list", {
  d1 <- make_boundary_block(c("A", "B"))
  d2 <- make_boundary_block(c("C", "D"), offset = 20L)

  expect_equal(bind_multidesign(d1)$x, d1$x)
  expect_equal(nrow(bind_multidesign(list(d1, d2))$x), 4L)
})

test_that("df_to_hyperdesign forwards the phase-1 contract", {
  input <- data.frame(
    subject = rep(c("one", "two"), each = 2L),
    entity = rep(c("A", "B"), 2L),
    condition = rep(c("x", "y"), 2L),
    v1 = 1:4,
    v2 = 5:8
  )

  hd <- df_to_hyperdesign(
    input,
    design_vars = c("entity", "condition"),
    x_vars = c("v1", "v2"),
    split_var = "subject",
    id = "entity",
    space = "common"
  )

  expect_equal(entity_id(hd), "entity")
  expect_equal(column_space(hd), "common")
  expect_equal(correspondence(hd)$global_ids, c("A", "B"))
})

test_that("block_indices remains independent of correspondence", {
  d1 <- make_boundary_block(c("A", "B"))
  d2 <- make_boundary_block(c("B", "C", "D"), offset = 20L)
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity")

  expect_equal(block_indices(hd), list(1:2, 3:4))
  expect_equal(block_indices(hd, byrow = TRUE), list(1:2, 3:5))
})
