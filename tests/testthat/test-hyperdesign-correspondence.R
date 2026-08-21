make_correspondence_block <- function(ids, ncol = 2L, id_name = "entity",
                                      column_design = NULL) {
  design <- data.frame(ids, stringsAsFactors = FALSE)
  names(design) <- id_name
  if (is.null(column_design)) {
    column_design <- data.frame(axis = seq_len(ncol))
  }
  multidesign(
    matrix(seq_len(max(1L, length(ids) * ncol)), nrow = length(ids), ncol = ncol),
    design,
    column_design
  )
}

test_that("unset hyperdesign contract preserves legacy behavior", {
  d1 <- make_correspondence_block(c("A", "B"))
  d2 <- make_correspondence_block(c("A", "C"))
  hd <- hyperdesign(list(d1, d2))

  expect_null(entity_id(hd))
  expect_null(column_space(hd))
  expect_false(has_correspondence(hd))
  expect_null(attr(hd, "entity_id", exact = TRUE))
  expect_null(attr(hd, "column_space", exact = TRUE))
  expect_null(attr(hd, "correspondence_assumption", exact = TRUE))
  expect_error(correspondence(hd), "no row-correspondence contract")
})

test_that("explicit correspondence represents partial overlap", {
  d1 <- make_correspondence_block(c("A", "B", "C"))
  d2 <- make_correspondence_block(c("B", "C", "D"))
  hd <- hyperdesign(
    list(first = d1, second = d2),
    id = "entity",
    space = "common"
  )
  cr <- correspondence(hd)

  expect_equal(entity_id(hd), "entity")
  expect_equal(column_space(hd), "common")
  expect_true(has_correspondence(hd))
  expect_s3_class(cr, "md_correspondence")
  expect_equal(cr$global_ids, c("A", "B", "C", "D"))
  expect_equal(cr$row_map$first, c(1L, 2L, 3L))
  expect_equal(cr$row_map$second, c(2L, 3L, 4L))
  expect_equal(cr$overlap$pair_n[1, 2], 2L)
  expect_true(cr$overlap$connected)

  for (nm in names(hd)) {
    expect_equal(cr$global_ids[cr$row_map[[nm]]], cr$ids[[nm]])
  }
  expect_equal(colSums(cr$overlap$incidence), cr$n_observed)
})

test_that("disconnected correspondence is reported rather than rejected", {
  hd <- hyperdesign(
    list(
      first = make_correspondence_block(c("A", "B")),
      second = make_correspondence_block(c("C", "D"))
    ),
    id = "entity"
  )

  expect_false(correspondence(hd)$overlap$connected)
})

test_that("global IDs use normalized first-seen order", {
  d1 <- make_correspondence_block(c(2L, 1L))
  d2 <- make_correspondence_block(c("1", "3"))
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity")

  expect_equal(correspondence(hd)$global_ids, c("2", "1", "3"))
  expect_type(correspondence(hd)$row_map[[1]], "integer")
  expect_type(d1$design$entity, "integer")
})

test_that("explicit entity IDs are validated at construction", {
  good <- make_correspondence_block(c("A", "B"))
  missing <- multidesign(matrix(1:4, 2), data.frame(other = 1:2))
  duplicate <- make_correspondence_block(c("A", "A"))
  with_na <- make_correspondence_block(c("A", NA_character_))
  list_id <- multidesign(
    matrix(1:4, 2),
    tibble::tibble(entity = list("A", "B"))
  )

  expect_error(hyperdesign(list(a = good, b = missing), id = "entity"), "missing")
  expect_error(hyperdesign(list(a = good, b = duplicate), id = "entity"), "duplicate")
  expect_error(hyperdesign(list(a = good, b = with_na), id = "entity"), "cannot contain NA")
  expect_error(hyperdesign(list(a = good, b = list_id), id = "entity"), "atomic vector")
  expect_error(hyperdesign(list(a = good), id = ".index"), "reserved")
  expect_error(hyperdesign(list(a = good), id = character()), "nonempty")
})

test_that("positional correspondence is explicit and requires equal rows", {
  d1 <- make_correspondence_block(c("A", "B"))
  d2 <- make_correspondence_block(c("X", "Y"))
  hd <- hyperdesign(list(one = d1, two = d2), positional = TRUE)
  cr <- correspondence(hd)

  expect_true(has_correspondence(hd))
  expect_null(entity_id(hd))
  expect_equal(cr$assumption, "positional")
  expect_equal(cr$global_ids, c("1", "2"))
  expect_equal(cr$row_map$one, 1:2)
  expect_error(
    hyperdesign(list(one = d1, two = make_correspondence_block(c("X"))), positional = TRUE),
    "equal row counts"
  )
  expect_error(
    hyperdesign(list(one = d1, two = d2), id = "entity", positional = TRUE),
    "mutually exclusive"
  )
  expect_error(hyperdesign(list(one = d1), positional = NA), "TRUE or FALSE")
})

test_that("column-space contracts validate dimensions and metadata", {
  cd <- data.frame(axis = c("x", "y"))
  d1 <- make_correspondence_block(c("A", "B"), column_design = cd)
  d2 <- make_correspondence_block(c("A", "B"), column_design = cd)
  d3 <- make_correspondence_block(c("A", "B"), ncol = 3L)
  d4 <- make_correspondence_block(
    c("A", "B"),
    column_design = data.frame(axis = c("u", "v"))
  )

  expect_equal(column_space(hyperdesign(list(a = d1, b = d2), space = "common")), "common")
  expect_equal(column_space(hyperdesign(list(a = d1, b = d3), space = "block")), "block")
  expect_error(hyperdesign(list(a = d1, b = d3), space = "common"), "equal column counts")
  expect_error(hyperdesign(list(a = d1, b = d4), space = "common"), "identical column designs")
  expect_error(hyperdesign(list(a = d1), space = "shared"), "common.*block")
})

test_that("contracted block names are unique and nonempty", {
  d1 <- make_correspondence_block(c("A", "B"))
  d2 <- make_correspondence_block(c("A", "C"))

  expect_error(
    hyperdesign(list(d1, d2), block_names = c("same", "same"), id = "entity"),
    "unique, nonempty"
  )
  expect_error(
    hyperdesign(stats::setNames(list(d1, d2), c("one", "")), id = "entity"),
    "unique, nonempty"
  )
})

test_that("correspondence invariants survive block permutation", {
  blocks <- list(
    one = make_correspondence_block(c("B", "A")),
    two = make_correspondence_block(c("A", "C")),
    three = make_correspondence_block(c("B", "C"))
  )
  cr1 <- correspondence(hyperdesign(blocks, id = "entity"))
  cr2 <- correspondence(hyperdesign(blocks[c(3, 1, 2)], id = "entity"))

  expect_equal(
    cr2$overlap$pair_n[names(blocks), names(blocks)],
    cr1$overlap$pair_n
  )
  expect_equal(sort(cr2$global_ids), sort(cr1$global_ids))
  expect_identical(cr2$overlap$connected, cr1$overlap$connected)
})

test_that("one-block correspondence is connected and revalidates current data", {
  block <- make_correspondence_block(c("A", "B", "C"))
  hd <- hyperdesign(list(only = block), id = "entity")
  cr <- correspondence(hd)

  expect_true(cr$overlap$connected)
  expect_equal(unname(cr$overlap$pair_n), matrix(3L, 1, 1))
  expect_equal(cr$global_ids[cr$row_map$only], c("A", "B", "C"))

  hd[[1]]$design$entity[[3]] <- "A"
  expect_error(correspondence(hd), "duplicate keys")
})

test_that("row order changes maps but not correspondence relationships", {
  blocks <- list(
    one = make_correspondence_block(c("A", "B", "C")),
    two = make_correspondence_block(c("B", "C", "D"))
  )
  original <- correspondence(hyperdesign(blocks, id = "entity"))
  permutation <- c(3L, 1L, 2L)
  blocks$one <- multidesign(
    blocks$one$x[permutation, , drop = FALSE],
    design(blocks$one)[permutation, , drop = FALSE],
    blocks$one$column_design
  )
  permuted <- correspondence(hyperdesign(blocks, id = "entity"))

  expect_equal(permuted$global_ids[permuted$row_map$one], c("C", "A", "B"))
  expect_equal(permuted$overlap$pair_n, original$overlap$pair_n)
  expect_equal(sort(permuted$global_ids), sort(original$global_ids))
})

test_that("NA values in x are not interpreted as correspondence masks", {
  d1 <- make_correspondence_block(c("A", "B"))
  d2 <- make_correspondence_block(c("A", "B"))
  d1$x[1, 1] <- NA_real_
  hd <- hyperdesign(list(one = d1, two = d2), id = "entity")

  expect_equal(correspondence(hd)$overlap$pair_n[1, 2], 2L)
  expect_true(correspondence(hd)$overlap$incidence["A", "one"])
})

test_that("print reports an opt-in contract without changing legacy output", {
  d1 <- make_correspondence_block(c("A", "B"))
  d2 <- make_correspondence_block(c("A", "C"))

  legacy_output <- crayon::strip_style(paste(
    testthat::capture_output(print(hyperdesign(list(d1, d2)))),
    collapse = "\n"
  ))
  expect_false(grepl("Contract:", legacy_output, fixed = TRUE))

  contracted <- hyperdesign(
    list(one = d1, two = d2),
    id = "entity",
    space = "common"
  )
  contract_output <- crayon::strip_style(paste(
    testthat::capture_output(print(contracted)),
    collapse = "\n"
  ))
  expect_match(contract_output, "Entity ID:\\s+entity")
  expect_match(contract_output, "Column space:\\s+common")
  expect_match(contract_output, "Global entities:\\s+3")
  expect_match(contract_output, "Overlap connected:\\s+TRUE")
})

test_that("optional sorting recomputes maps without changing first-seen defaults", {
  blocks <- list(
    one = make_correspondence_block(c("B", "A")),
    two = make_correspondence_block(c("A", "C")),
    three = make_correspondence_block(c("C", "B"))
  )
  hd <- hyperdesign(blocks, id = "entity")
  before <- hd

  first_seen <- correspondence(hd)
  sorted <- correspondence(hd, sort_ids = TRUE)

  expect_equal(first_seen$global_ids, c("B", "A", "C"))
  expect_equal(sorted$global_ids, c("A", "B", "C"))
  for (nm in names(blocks)) {
    expect_equal(sorted$global_ids[sorted$row_map[[nm]]], sorted$ids[[nm]])
  }
  expect_equal(sorted$overlap$pair_n, first_seen$overlap$pair_n)
  expect_identical(sorted$overlap$connected, first_seen$overlap$connected)
  expect_identical(hd, before)
  expect_equal(correspondence(hd)$global_ids, first_seen$global_ids)
  expect_error(correspondence(hd, sort_ids = NA), "TRUE or FALSE")
})

test_that("sorted correspondence is stable under block permutations", {
  blocks <- list(
    one = make_correspondence_block(c("D", "A", "B")),
    two = make_correspondence_block(c("C", "A")),
    three = make_correspondence_block(c("B", "E", "C"))
  )
  reference <- correspondence(hyperdesign(blocks, id = "entity"), sort_ids = TRUE)

  set.seed(20260820)
  for (iteration in seq_len(20L)) {
    permutation <- sample(seq_along(blocks))
    candidate <- correspondence(
      hyperdesign(blocks[permutation], id = "entity"),
      sort_ids = TRUE
    )
    expect_equal(candidate$global_ids, reference$global_ids)
    expect_equal(
      candidate$overlap$pair_n[names(blocks), names(blocks)],
      reference$overlap$pair_n
    )
    for (nm in names(blocks)) {
      expect_equal(candidate$row_map[[nm]], reference$row_map[[nm]])
    }
  }
})

test_that("sorting supports positional and normalized mixed-type IDs", {
  positional <- hyperdesign(
    list(
      one = make_correspondence_block(c("x", "y", "z")),
      two = make_correspondence_block(c("a", "b", "c"))
    ),
    positional = TRUE
  )
  expect_equal(
    correspondence(positional, sort_ids = TRUE)$global_ids,
    c("1", "2", "3")
  )

  mixed <- hyperdesign(
    list(
      one = make_correspondence_block(c(10L, 2L)),
      two = make_correspondence_block(c("2", "1"))
    ),
    id = "entity"
  )
  sorted <- correspondence(mixed, sort_ids = TRUE)
  expect_equal(sorted$global_ids, c("1", "10", "2"))
  expect_equal(sorted$global_ids[sorted$row_map$one], c("10", "2"))
})
