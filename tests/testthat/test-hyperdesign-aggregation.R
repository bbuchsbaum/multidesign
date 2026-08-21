make_duplicate_entity_block <- function(with_cells = FALSE,
                                        conflicting_metadata = FALSE) {
  x <- rbind(
    c(1, 100, NA, 40),
    c(3, 5, 7, 50),
    c(9, 10, 11, 12)
  )
  design <- data.frame(
    entity = c("A", "A", "B"),
    curve = c("jaw", if (conflicting_metadata) "orbit" else "jaw", "orbit"),
    replicate = c(1L, 1L, 2L)
  )
  cells <- if (with_cells) {
    rbind(
      c(TRUE, FALSE, TRUE, FALSE),
      c(TRUE, TRUE, FALSE, FALSE),
      c(FALSE, FALSE, FALSE, FALSE)
    )
  } else {
    NULL
  }
  multidesign(
    x,
    design,
    data.frame(axis = paste0("v", seq_len(4))),
    cells = cells
  )
}

test_that("duplicate entity aggregation remains opt-in", {
  block <- make_duplicate_entity_block()

  expect_error(
    hyperdesign(list(one = block), id = "entity"),
    "duplicate keys"
  )
  expect_error(
    hyperdesign(list(one = block), aggregate = "mean"),
    "requires an explicit entity"
  )
  expect_error(
    hyperdesign(list(one = block), id = "entity", aggregate = "median"),
    "mean.*function"
  )

  aggregated <- hyperdesign(
    list(one = block),
    id = "entity",
    space = "common",
    aggregate = "mean"
  )
  expect_equal(design(aggregated$one)$entity, c("A", "B"))
  expect_equal(aggregated$one$x[1, ], c(2, 52.5, NA, 45))
  expect_equal(aggregated$one$x[2, ], block$x[3, ])
  expect_equal(correspondence(aggregated)$global_ids, c("A", "B"))
})

test_that("custom duplicate aggregation is scalar and preserves singleton rows", {
  block <- make_duplicate_entity_block()
  aggregated <- hyperdesign(
    list(specimen = block),
    id = "entity",
    aggregate = function(values) max(values, na.rm = TRUE)
  )

  expect_equal(aggregated$specimen$x[1, ], c(3, 100, 7, 50))
  expect_equal(aggregated$specimen$x[2, ], block$x[3, ])

  expect_error(
    hyperdesign(
      list(specimen = block),
      id = "entity",
      aggregate = function(values) c(min(values), max(values))
    ),
    "block specimen.*entity.*A.*coordinate 1.*one numeric"
  )
})

test_that("duplicate aggregation rejects ambiguous non-key metadata", {
  block <- make_duplicate_entity_block(conflicting_metadata = TRUE)

  expect_error(
    hyperdesign(list(one = block), id = "entity", aggregate = "mean"),
    "entity `A`.*block `one`.*curve"
  )
})

test_that("duplicate aggregation is mask-aware without treating NA as a mask", {
  block <- make_duplicate_entity_block(with_cells = TRUE)
  aggregated <- hyperdesign(
    list(one = block),
    id = "entity",
    space = "common",
    aggregate = "mean"
  )$one

  expect_equal(aggregated$x[1, 1:2], c(2, 5))
  expect_true(is.na(aggregated$x[1, 3]))
  expect_true(cell_mask(aggregated)[1, 3])
  expect_true(is.na(aggregated$x[1, 4]))
  expect_false(cell_mask(aggregated)[1, 4])

  expect_equal(aggregated$x[2, ], block$x[3, ])
  expect_false(any(cell_mask(aggregated)[2, ]))
})

test_that("df_to_hyperdesign forwards explicit duplicate aggregation", {
  data <- data.frame(
    specimen = c("one", "one", "one"),
    entity = c("A", "A", "B"),
    curve = c("jaw", "jaw", "orbit"),
    x = c(1, 3, 9),
    y = c(2, 4, 10)
  )
  hd <- df_to_hyperdesign(
    data,
    design_vars = c("entity", "curve"),
    x_vars = c("x", "y"),
    split_var = "specimen",
    id = "entity",
    space = "common",
    aggregate = "mean"
  )

  expect_equal(hd$one$x, rbind(c(2, 3), c(9, 10)), ignore_attr = TRUE)
  expect_equal(correspondence(hd)$global_ids, c("A", "B"))
})

test_that("hyperdesign summaries preserve only truthful correspondence", {
  column_design <- data.frame(axis = c("x", "y"))
  one <- multidesign(
    matrix(1:6, nrow = 3),
    data.frame(entity = c("A", "B", "C"), curve = c("jaw", "jaw", "orbit")),
    column_design
  )
  two <- multidesign(
    matrix(7:12, nrow = 3),
    data.frame(entity = c("A", "C", "D"), curve = c("jaw", "orbit", "jaw")),
    column_design
  )
  hd <- hyperdesign(list(one = one, two = two), id = "entity", space = "common")

  by_entity <- summarize_by(hd, entity)
  expect_equal(entity_id(by_entity), "entity")
  expect_equal(correspondence(by_entity)$global_ids, c("A", "B", "C", "D"))

  by_curve <- summarize_by(hd, curve)
  expect_null(entity_id(by_curve))
  expect_false(has_correspondence(by_curve))
  expect_equal(column_space(by_curve), "common")

  positional <- hyperdesign(list(one = one, two = one), positional = TRUE, space = "common")
  positional_summary <- summarize_by(positional, curve)
  expect_false(has_correspondence(positional_summary))
  expect_equal(column_space(positional_summary), "common")
})
