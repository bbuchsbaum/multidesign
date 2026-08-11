.adversarial_frame <- function(seed, lazy_block = FALSE) {
  set.seed(seed)
  n <- sample(12:24, 1L)
  fac1 <- factor(sample(c("A", "B", "C"), n, replace = TRUE), levels = c("A", "B", "C"))
  fac2 <- factor(sample(c("old", "new"), n, replace = TRUE), levels = c("old", "new"))
  x <- rnorm(n)
  block_values <- matrix(rnorm(n * 3L), nrow = n)
  block_data <- if (lazy_block) {
    fmridataset::counting_source(fmridataset::memory_source(block_values))
  } else {
    block_values
  }
  frame <- fmridataset::fmri_frame(
    assays = list(beta = matrix(rnorm(n * 5L), nrow = n)),
    observations = fmridataset::axis_frame(
      data.frame(
        .obs_id = paste0("seed-", seed, "-", seq_len(n)),
        subject_id = factor(sample(paste0("s", 1:5), n, replace = TRUE)),
        Fac1 = fac1,
        Fac2 = fac2,
        x = x
      ),
      blocks = list(
        features = fmridataset::axis_block(
          block_data,
          components = data.frame(.component_id = c("u", "v", "w"))
        )
      )
    )
  )
  list(frame = frame, block = block_values)
}

test_that("random ragged designs match dense reference calculations", {
  for (seed in 1:20) {
    fixture <- .adversarial_frame(seed)
    compiled <- compile_design(
      fixture$frame,
      design_spec(~ Fac1 * Fac2 + x + mv(features, c("w", "u")))
    )
    data <- fmridataset::observations(fixture$frame)
    data$mv__features__w <- fixture$block[, 3L]
    data$mv__features__u <- fixture$block[, 1L]
    reference <- stats::model.matrix(
      ~ Fac1 * Fac2 + x + mv__features__w + mv__features__u,
      data = data
    )
    expect_equal(model_matrix(compiled), reference, tolerance = 0, info = seed)

    order <- sample(seq_len(nrow(fixture$frame)))
    reordered <- compile_design(fixture$frame[order, c(5L, 1L)], compiled$spec)
    expect_equal(
      unname(model_matrix(reordered)),
      unname(reference[order, , drop = FALSE]),
      tolerance = 0,
      ignore_attr = TRUE,
      info = paste("permutation", seed)
    )
  }
})

test_that("lazy and in-memory block backends are design-equivalent", {
  memory <- .adversarial_frame(91L, lazy_block = FALSE)
  lazy <- .adversarial_frame(91L, lazy_block = TRUE)
  spec <- design_spec(~ Fac1 + mv(features, c("v", "u")))
  memory_design <- compile_design(memory$frame, spec)
  lazy_design <- compile_design(lazy$frame, spec)

  expect_equal(model_matrix(lazy_design), model_matrix(memory_design), tolerance = 0)
  expect_identical(design_digest(lazy_design), design_digest(memory_design))
  source <- fmridataset::axis_block_data(fmridataset::obs_blocks(lazy$frame)$features)
  expect_equal(fmridataset::source_counts(source)$reads, 1)
})

test_that("hostile component IDs and overlapping mv terms remain unambiguous", {
  frame <- .adversarial_frame(44L)$frame
  values <- matrix(seq_len(nrow(frame) * 3L), nrow = nrow(frame))
  frame$observations$blocks$hostile <- fmridataset::axis_block(
    values,
    components = data.frame(.component_id = c("a-b", "a b", "x:y"))
  )
  compiled <- compile_design(
    frame,
    design_spec(
      ~ 0 + mv(hostile, "a b") + mv(hostile, c("a-b", "a b", "x:y"))
    )
  )

  expect_identical(ncol(model_matrix(compiled)), 3L)
  expect_identical(component_data(compiled)$component_id, c("a b", "a-b", "x:y"))
  expect_identical(nrow(compiled$component_map), 3L)
  expect_equal(
    unname(model_matrix(compiled)),
    values[, c(2L, 1L, 3L), drop = FALSE],
    tolerance = 0,
    ignore_attr = TRUE
  )
  restored <- apply_design(compiled, frame[c(3L, 1L, 2L), ])
  expect_equal(
    unname(model_matrix(restored)),
    values[c(3L, 1L, 2L), c(2L, 1L, 3L), drop = FALSE],
    tolerance = 0,
    ignore_attr = TRUE
  )
})

test_that("adversarial selectors and cache bounds fail deterministically", {
  frame <- .adversarial_frame(8L)$frame
  chosen <- c("u", "missing")
  expect_error(
    design_input_digest(frame, design_spec(~ mv(features, chosen))),
    "selection is invalid"
  )
  expect_error(design_cache(0), "positive integer")
  expect_error(design_cache(1.5), "positive integer")
  expect_error(clear_design_cache(new.env()), "design_cache")
})
