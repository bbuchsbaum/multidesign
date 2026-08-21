make_contract_fold_block <- function(ids, run = seq_along(ids), offset = 0) {
  multidesign(
    matrix(seq_len(length(ids) * 2L) + offset, nrow = length(ids)),
    data.frame(entity = ids, run = run),
    data.frame(axis = c("x", "y"))
  )
}

test_that("fold_over on the entity key holds out global entities", {
  d1 <- make_contract_fold_block(c("A", "B", "C"))
  d2 <- make_contract_fold_block(c("A", "B", "D"), offset = 20)
  hd <- hyperdesign(
    list(one = d1, two = d2),
    id = "entity",
    space = "common"
  )

  folds <- fold_over(hd, entity)
  expect_length(folds, 4L)

  held <- vapply(folds, function(fold) fold$held_out$entity, character(1))
  expect_equal(held, c("A", "B", "C", "D"))
  materialized <- lapply(seq_along(folds), function(i) folds[[i]])
  for (fold in materialized) {
    entity <- fold$held_out$entity
    expect_s3_class(fold$analysis, "hyperdesign")
    expect_s3_class(fold$assessment, "hyperdesign")
    expect_equal(entity_id(fold$analysis), "entity")
    expect_equal(entity_id(fold$assessment), "entity")
    expect_false(any(vapply(
      fold$analysis,
      function(block) entity %in% as.character(block$design$entity),
      logical(1)
    )))
    expect_true(all(vapply(
      fold$assessment,
      function(block) identical(as.character(block$design$entity), entity),
      logical(1)
    )))
  }
  expect_output(print(folds), "Assessment Set")
})

test_that("non-entity folds remain block-local and preserve explicit contract", {
  d1 <- make_contract_fold_block(c("A", "B", "C", "D"), c(1, 1, 2, 2))
  d2 <- make_contract_fold_block(c("A", "B", "C", "D"), c(1, 1, 2, 2), 20)
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity", space = "common")

  folds <- fold_over(hd, run)
  expect_length(folds, 4L)
  first <- folds[[1]]
  expect_s3_class(first$assessment, "multidesign")
  expect_equal(entity_id(first$analysis), "entity")
  expect_true("entity" %in% names(first$assessment$design))
  expect_equal(nrow(first$analysis[[1]]$x), 2L)
  expect_equal(nrow(first$analysis[[2]]$x), 4L)
})

test_that("entity folds drop an emptied analysis block but retain the problem", {
  d1 <- make_contract_fold_block("A")
  d2 <- make_contract_fold_block(c("A", "B"), offset = 20)
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity", space = "common")

  fold_a <- fold_over(hd, entity)[[1]]
  expect_equal(names(fold_a$analysis), "two")
  expect_equal(as.character(fold_a$analysis[[1]]$design$entity), "B")
  expect_equal(names(fold_a$assessment), c("one", "two"))
})

test_that("cv_rows preserves explicit contracts on analysis and assessment", {
  d1 <- make_contract_fold_block(c("A", "B", "C"))
  d2 <- make_contract_fold_block(c("A", "B", "D"), offset = 20)
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity", space = "common")

  fold <- cv_rows(
    hd,
    rows = list(list(one = 1L, two = 2L)),
    preserve_row_ids = TRUE
  )[[1]]

  expect_equal(entity_id(fold$analysis), "entity")
  expect_equal(entity_id(fold$assessment), "entity")
  expect_equal(column_space(fold$analysis), "common")
  expect_equal(column_space(fold$assessment), "common")
  expect_equal(fold$assessment[[1]]$design$entity, "A")
  expect_equal(fold$assessment[[2]]$design$entity, "B")
  expect_equal(fold$held_out$row_ids$one, 1L)
  expect_equal(fold$held_out$row_ids$two, 2L)
})

test_that("positional cv_rows accepts synchronized and rejects divergent rows", {
  d1 <- make_contract_fold_block(c("A", "B", "C"))
  d2 <- make_contract_fold_block(c("X", "Y", "Z"), offset = 20)
  hd <- hyperdesign(list(one = d1, two = d2), positional = TRUE, space = "common")

  synchronized <- cv_rows(
    hd,
    rows = list(list(one = 2L, two = 2L))
  )[[1]]
  expect_true(has_correspondence(synchronized$analysis))
  expect_true(has_correspondence(synchronized$assessment))

  divergent <- cv_rows(
    hd,
    rows = list(list(one = 1L, two = 2L))
  )
  expect_error(divergent[[1]], "identical retained row positions")

  block_local <- fold_over(hd, run)
  expect_error(block_local[[1]], "identical retained row positions")
})

test_that("entity ID columns survive folds without preserve_row_ids", {
  d1 <- make_contract_fold_block(c("A", "B"))
  d2 <- make_contract_fold_block(c("A", "B"), offset = 20)
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity")

  fold <- fold_over(hd, entity)[[1]]
  expect_true(all(vapply(
    fold$assessment,
    function(block) "entity" %in% names(block$design),
    logical(1)
  )))
  expect_false(any(vapply(
    fold$assessment,
    function(block) ".orig_index" %in% names(block$design),
    logical(1)
  )))
})

test_that("entity folds reject a fold with no analysis data", {
  d1 <- make_contract_fold_block("A")
  d2 <- make_contract_fold_block("A", offset = 20)
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity")

  folds <- fold_over(hd, entity)
  expect_error(folds[[1]], "at least one analysis block")
})
