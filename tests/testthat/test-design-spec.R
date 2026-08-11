.design_frame_fixture <- function() {
  observations <- data.frame(
    .obs_id = paste0("obs-", 1:8),
    subject_id = factor(rep(c("s01", "s02", "s03", "s04"), each = 2L)),
    stimulus_id = c("i1", "i2", "i2", "i3", "i1", "i3", "i3", "i1"),
    Fac1 = factor(rep(c("A", "B"), 4L), levels = c("A", "B")),
    Fac2 = factor(c("old", "new", "new", "old", "old", "new", "old", "new")),
    age = c(61, 61, 67, 67, 72, 72, 64, 64)
  )
  motion <- fmridataset::axis_block(
    matrix(seq_len(16) / 10, nrow = 8L),
    components = data.frame(.component_id = c("translation", "rotation")),
    role = "confound"
  )
  stimulus <- list(
    key = "stimulus_id",
    data = data.frame(
      stimulus_id = c("i1", "i2", "i3"),
      category = c("face", "scene", "object")
    ),
    blocks = list(
      visual_pca = fmridataset::axis_block(
        matrix(c(.1, .2, .3, .4, .5, .6, .7, .8, .9), nrow = 3L),
        components = data.frame(.component_id = c("PC01", "PC02", "PC03"))
      )
    )
  )
  fmridataset::fmri_frame(
    assays = list(beta = matrix(seq_len(40), nrow = 8L)),
    observations = fmridataset::axis_frame(
      observations,
      blocks = list(motion = motion)
    ),
    entities = list(stimulus = stimulus)
  )
}

test_that("design_spec preserves formulas and validates contracts", {
  spec <- design_spec(
    fixed = ~ Fac1 + mv(motion, 1:2),
    random = ~ 1 | subject_id
  )
  expect_s3_class(spec, "design_spec")
  expect_identical(spec$na_action, "fail")
  expect_identical(environment(spec$fixed), environment())
  expect_error(design_spec(1), "one-sided")
  expect_error(design_spec(y ~ Fac1), "one-sided")
  expect_error(design_spec(~Fac1, random = y ~ subject_id), "one-sided")
  expect_error(mv(motion), "formula special")
})

test_that("compile_design matches an explicit dense reference", {
  frame <- .design_frame_fixture()
  spec <- design_spec(
    fixed = ~ Fac1 * Fac2 + age + mv(motion) + mv(stimulus.visual_pca, 1:2),
    random = ~ 1 + Fac1 | subject_id
  )
  compiled <- compile_design(frame, spec)
  observations <- fmridataset::observations(frame)
  motion <- fmridataset::axis_block_data(fmridataset::obs_blocks(frame)$motion)
  stimuli <- frame$entities$stimulus
  stimulus_index <- match(observations$stimulus_id, stimuli$data$stimulus_id)
  visual <- fmridataset::axis_block_data(stimuli$blocks$visual_pca)[stimulus_index, 1:2]
  reference_data <- cbind(
    observations,
    mv__motion__translation = motion[, 1L],
    mv__motion__rotation = motion[, 2L],
    mv__stimulus.visual_pca__PC01 = visual[, 1L],
    mv__stimulus.visual_pca__PC02 = visual[, 2L]
  )
  reference <- stats::model.matrix(
    ~ Fac1 * Fac2 + age +
      mv__motion__translation + mv__motion__rotation +
      mv__stimulus.visual_pca__PC01 + mv__stimulus.visual_pca__PC02,
    data = reference_data
  )

  expect_equal(model_matrix(compiled), reference)
  expect_identical(compiled$observation_ids, fmridataset::observation_ids(frame))
  expect_identical(names(grouping_data(compiled)), "subject_id")
  expect_identical(grouping_data(compiled)$subject_id, observations$subject_id)
  expect_setequal(
    na.omit(unique(term_data(compiled)$source_block)),
    c("motion", "stimulus.visual_pca")
  )
  expect_setequal(
    na.omit(unique(term_data(compiled)$component_id)),
    c("translation", "rotation", "PC01", "PC02")
  )
})

test_that("entity blocks lift lazily and track stable component IDs", {
  frame <- .design_frame_fixture()
  compiled <- compile_design(
    frame,
    design_spec(~ 0 + mv(stimulus.visual_pca, c("PC03", "PC01")))
  )
  entity <- frame$entities$stimulus
  index <- match(
    fmridataset::observations(frame)$stimulus_id,
    entity$data$stimulus_id
  )
  expected <- fmridataset::axis_block_data(entity$blocks$visual_pca)[index, c(3L, 1L)]
  colnames(expected) <- colnames(model_matrix(compiled))

  expect_equal(model_matrix(compiled), expected, ignore_attr = TRUE)
  expect_identical(
    term_data(compiled)$component_id,
    c("PC03", "PC01")
  )
  expect_true(all(compiled$component_map$alignment == "stimulus"))
})

test_that("compiled designs remain correct after synchronized row views", {
  frame <- .design_frame_fixture()
  view <- fmridataset::filter_obs(frame, Fac1 == "A")
  compiled <- compile_design(
    view,
    design_spec(~ Fac2 + mv(motion, "rotation"))
  )
  data <- fmridataset::observations(view)
  rotation <- fmridataset::axis_block_data(fmridataset::obs_blocks(view)$motion)[, 2L]
  reference <- stats::model.matrix(
    ~ Fac2 + mv__motion__rotation,
    data = transform(data, mv__motion__rotation = rotation)
  )
  expect_equal(model_matrix(compiled), reference)
  expect_identical(compiled$observation_ids, fmridataset::observation_ids(view))
})

test_that("compile_design rejects unresolved blocks keys and components", {
  frame <- .design_frame_fixture()
  expect_error(
    compile_design(frame, design_spec(~ mv(unknown))),
    "Unknown multivariate block"
  )
  expect_error(
    compile_design(frame, design_spec(~ mv(motion, 99))),
    "selection is invalid"
  )
  broken <- frame
  broken$entities$stimulus$data$stimulus_id[1L] <- "missing"
  expect_error(
    compile_design(broken, design_spec(~ mv(stimulus.visual_pca))),
    "unresolved"
  )
})

test_that("na policies retain explicit observation identity", {
  frame <- .design_frame_fixture()
  frame$observations$data$age[3L] <- NA_real_
  expect_error(
    compile_design(frame, design_spec(~age, na_action = "fail")),
    "missing values"
  )
  compiled <- compile_design(frame, design_spec(~age, na_action = "omit"))
  expect_identical(
    compiled$observation_ids,
    fmridataset::observation_ids(frame)[-3L]
  )
})
