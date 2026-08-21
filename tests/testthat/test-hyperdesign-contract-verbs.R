make_contract_verb_block <- function(ids, keep, offset = 0, column_design = NULL) {
  if (is.null(column_design)) {
    column_design <- data.frame(
      axis = paste0("v", 1:3),
      region = c("left", "left", "right")
    )
  }
  multidesign(
    matrix(seq_len(length(ids) * nrow(column_design)) + offset,
           nrow = length(ids)),
    data.frame(entity = ids, keep = keep),
    column_design
  )
}

test_that("subset rebuilds explicit correspondence from filtered blocks", {
  d1 <- make_contract_verb_block(c("A", "B", "C"), c(TRUE, TRUE, FALSE))
  d2 <- make_contract_verb_block(c("B", "C", "D"), c(FALSE, TRUE, TRUE), 20)
  hd <- hyperdesign(
    list(one = d1, two = d2),
    id = "entity",
    space = "common"
  )

  filtered <- subset(hd, keep)
  cr <- correspondence(filtered)

  expect_equal(entity_id(filtered), "entity")
  expect_equal(column_space(filtered), "common")
  expect_equal(cr$global_ids, c("A", "B", "C", "D"))
  expect_false(cr$overlap$connected)
  expect_equal(attr(filtered, "hdes")$nr, c(2L, 2L))
  expect_equal(attr(filtered, "common_vars"), c("entity", "keep"))
})

test_that("subset preserves only synchronized positional correspondence", {
  d1 <- make_contract_verb_block(c("A", "B", "C"), c(TRUE, FALSE, TRUE))
  d2 <- make_contract_verb_block(c("X", "Y", "Z"), c(TRUE, FALSE, TRUE), 20)
  hd <- hyperdesign(list(one = d1, two = d2), positional = TRUE)

  filtered <- subset(hd, keep)
  expect_equal(correspondence(filtered)$global_ids, c("1", "2"))
  expect_true(has_correspondence(filtered))

  d2_bad <- make_contract_verb_block(c("X", "Y", "Z"), c(TRUE, TRUE, FALSE), 20)
  hd_bad <- hyperdesign(list(one = d1, two = d2_bad), positional = TRUE)
  expect_error(
    subset(hd_bad, keep),
    "identical retained row positions"
  )
})

test_that("select_variables preserves and revalidates common space", {
  d1 <- make_contract_verb_block(c("A", "B"), c(TRUE, TRUE))
  d2 <- make_contract_verb_block(c("A", "C"), c(TRUE, TRUE), 20)
  hd <- hyperdesign(
    list(one = d1, two = d2),
    id = "entity",
    space = "common"
  )

  selected <- select_variables(hd, region == "left")
  expect_equal(
    unname(vapply(selected, function(block) ncol(block$x), integer(1))),
    c(2L, 2L)
  )
  expect_equal(entity_id(selected), "entity")
  expect_equal(column_space(selected), "common")
  expect_equal(correspondence(selected)$global_ids, c("A", "B", "C"))

  hd[[2]]$column_design$region <- c("left", "right", "right")
  expect_error(
    select_variables(hd, region == "left"),
    "equal column counts|identical column designs"
  )
})

test_that("init_transform preserves names, metadata, and contract", {
  d1 <- make_contract_verb_block(c("A", "B", "C"), rep(TRUE, 3))
  d2 <- make_contract_verb_block(c("A", "C", "D"), rep(TRUE, 3), 20)
  hd <- hyperdesign(
    list(one = d1, two = d2),
    id = "entity",
    space = "common"
  )
  transformed <- init_transform(hd, multivarious::center())

  expect_equal(names(transformed), c("one", "two"))
  expect_equal(entity_id(transformed), "entity")
  expect_equal(column_space(transformed), "common")
  expect_equal(transformed[[1]]$column_design, d1$column_design)
  expect_equal(correspondence(transformed)$global_ids, c("A", "B", "C", "D"))
  expect_named(attr(transformed, "preproc"), c("one", "two"))
})
