.cache_frame_fixture <- function() {
  block_source <- fmridataset::counting_source(fmridataset::memory_source(
    matrix(c(.1, .2, .3, .4, .5, .6), nrow = 3L, byrow = TRUE)
  ))
  assay_source <- fmridataset::counting_source(fmridataset::memory_source(
    matrix(seq_len(18L), nrow = 6L)
  ))
  stimulus <- fmridataset::entity_frame(
    data = data.frame(stimulus_id = c("i1", "i2", "i3")),
    key = "stimulus_id",
    blocks = list(
      embedding = fmridataset::axis_block(
        block_source,
        components = data.frame(.component_id = c("e1", "e2"))
      )
    )
  )
  frame <- fmridataset::fmri_frame(
    assays = list(beta = assay_source),
    observations = data.frame(
      .obs_id = paste0("obs-", seq_len(6L)),
      subject_id = factor(rep(c("s1", "s2", "s3"), each = 2L)),
      stimulus_id = rep(c("i1", "i2", "i3"), each = 2L),
      condition = factor(rep(c("A", "B"), 3L)),
      x = seq_len(6L)
    ),
    entities = list(stimulus = stimulus),
    relations = list(
      observation_stimulus = fmridataset::key_relation(
        "stimulus_id", target = "stimulus"
      )
    )
  )
  list(
    frame = frame,
    block_source = block_source,
    assay_source = assay_source
  )
}

.cache_design_spec <- function() {
  design_spec(
    ~ condition + scale(x) + mv(stimulus.embedding, "e1"),
    ~ 1 | subject_id
  )
}

test_that("design digests track semantic inputs but not imaging layout", {
  fixture <- .cache_frame_fixture()
  spec <- .cache_design_spec()
  input <- design_input_digest(fixture$frame, spec)
  compiled <- compile_design(fixture$frame, spec)

  expect_match(input, "^[0-9a-f]{64}$")
  expect_match(design_digest(compiled), "^[0-9a-f]{64}$")
  expect_identical(
    design_digest(unserialize(serialize(compiled, NULL))),
    design_digest(compiled)
  )
  expect_identical(
    design_input_digest(fixture$frame[, 1L], spec),
    input
  )
  expect_equal(fmridataset::source_counts(fixture$assay_source)$bytes, 0)

  changed <- fixture$frame
  changed$observations$data$x[[1L]] <- 99
  expect_false(identical(design_input_digest(changed, spec), input))
  expect_false(identical(
    design_digest(compile_design(changed, spec)),
    design_digest(compiled)
  ))
  reordered <- fixture$frame[c(6L, 1:5), ]
  expect_false(identical(design_input_digest(reordered, spec), input))

  unused <- fixture$frame
  unused$observations$data$unused <- seq_len(nrow(unused))
  expect_identical(design_input_digest(unused, spec), input)
  chosen <- "e1"
  dynamic <- design_spec(~ mv(stimulus.embedding, chosen))
  dynamic_e1 <- design_input_digest(fixture$frame, dynamic)
  chosen <- "e2"
  expect_false(identical(
    design_input_digest(fixture$frame, dynamic),
    dynamic_e1
  ))
})

test_that("design caches are bounded isolated and dependency aware", {
  fixture <- .cache_frame_fixture()
  spec <- .cache_design_spec()
  cache <- design_cache(max_entries = 2L)
  first <- compile_design(fixture$frame, spec, cache = cache)
  first_counts <- fmridataset::source_counts(fixture$block_source)
  first$model_matrix[[1L]] <- 999
  second <- compile_design(fixture$frame[, 2:3], spec, cache = cache)

  expect_false(identical(second$model_matrix[[1L]], 999))
  expect_identical(
    fmridataset::source_counts(fixture$block_source)$reads,
    first_counts$reads
  )
  expect_identical(
    design_cache_info(cache)[c("entries", "hits", "misses")],
    list(entries = 1L, hits = 1L, misses = 1L)
  )
  expect_equal(fmridataset::source_counts(fixture$assay_source)$bytes, 0)

  changed <- fixture$frame
  changed$observations$data$x[[1L]] <- 100
  compile_design(changed, spec, cache = cache)
  expect_identical(design_cache_info(cache)$misses, 2L)
  clear_design_cache(cache)
  expect_identical(design_cache_info(cache)$entries, 0L)
})

test_that("compiled frame folds fit transforms on analysis observations only", {
  fixture <- .cache_frame_fixture()
  folds <- compile_design_folds(
    fixture$frame,
    .cache_design_spec(),
    assessment = list(s1 = 1:2, s2 = 3:4)
  )

  expect_s3_class(folds, "compiled_design_folds")
  expect_length(folds, 2L)
  first <- folds[[1L]]
  expect_identical(first$analysis_design$observation_ids, paste0("obs-", 3:6))
  expect_identical(first$assessment_design$observation_ids, paste0("obs-", 1:2))
  expected_scaled <- (1:2 - mean(3:6)) / stats::sd(3:6)
  expect_equal(
    unname(model_matrix(first$assessment_design)[, "scale(x)"]),
    expected_scaled,
    tolerance = 0
  )
  expect_false(isTRUE(all.equal(
    expected_scaled,
    as.numeric(scale(1:2))
  )))
  index <- design_fold_data(folds)
  expect_identical(names(index), c(".fold", ".obs_id", "role"))
  expect_equal(as.list(table(index$role)), list(analysis = 8L, assessment = 4L))
  expect_equal(fmridataset::source_counts(fixture$assay_source)$bytes, 0)
})

test_that("compiled frame folds validate assessment selectors", {
  frame <- .cache_frame_fixture()$frame
  spec <- .cache_design_spec()
  expect_error(compile_design_folds(frame, spec, list(integer())), "non-empty")
  expect_error(compile_design_folds(frame, spec, list(c(1L, 1L))), "unique")
  expect_error(compile_design_folds(frame, spec, list(1:6)), "analysis")
  expect_error(compile_design_folds(frame, spec, list("unknown")), "unknown")
})
