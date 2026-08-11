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
  stimulus <- fmridataset::entity_frame(
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
  subject <- fmridataset::entity_frame(
    key = "subject_id",
    data = data.frame(
      subject_id = c("s01", "s02", "s03", "s04"),
      age = c(61, 67, 72, 64),
      site = factor(c("north", "south", "north", "south"))
    )
  )
  fmridataset::fmri_frame(
    assays = list(beta = matrix(seq_len(40), nrow = 8L)),
    observations = fmridataset::axis_frame(
      observations,
      blocks = list(motion = motion)
    ),
    entities = list(stimulus = stimulus, subject = subject),
    relations = list(
      observation_stimulus = fmridataset::key_relation(
        "stimulus_id", target = "stimulus"
      ),
      observation_subject = fmridataset::key_relation(
        "subject_id", target = "subject"
      )
    )
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
    "missing values"
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

test_that("resolved entity scalars participate in fixed and grouping formulas", {
  frame <- .design_frame_fixture()
  compiled <- compile_design(
    frame,
    design_spec(
      ~ Fac1 + subject.age + subject.site,
      random = ~ 1 | subject.subject_id
    )
  )
  resolved <- fmridataset::observations(frame, resolve = TRUE)
  reference <- stats::model.matrix(
    ~ Fac1 + subject.age + subject.site,
    data = resolved
  )

  expect_equal(model_matrix(compiled), reference)
  expect_identical(
    grouping_data(compiled)$subject.subject_id,
    resolved$subject.subject_id
  )
})

test_that("entity mv blocks use relation-aware lazy lifting after row selection", {
  frame <- .design_frame_fixture()
  stimulus <- fmridataset::entity(frame, "stimulus")
  values <- fmridataset::axis_block_data(
    fmridataset::entity_blocks(stimulus)$visual_pca
  )
  counted <- fmridataset::counting_source(fmridataset::memory_source(
    values
  ))
  stimulus$blocks$visual_pca$data <- counted
  frame$entities$stimulus <- stimulus
  view <- frame[c(1L, 5L, 8L), ]

  compiled <- compile_design(
    view,
    design_spec(~ 0 + mv(stimulus.visual_pca, c("PC03", "PC01")))
  )

  expect_equal(
    unname(model_matrix(compiled)),
    values[rep(1L, 3L), c(3L, 1L), drop = FALSE],
    ignore_attr = TRUE
  )
  counts <- fmridataset::source_counts(counted)
  expect_identical(counts$reads, 1)
  expect_identical(counts$bytes, 2 * 8)
})

test_that("term attribution uses exact generated variables rather than prefixes", {
  frame <- .design_frame_fixture()
  block <- fmridataset::axis_block(
    matrix(seq_len(16L), nrow = 8L),
    components = data.frame(.component_id = c("x", "x1"))
  )
  frame$observations$blocks$prefix <- block
  compiled <- compile_design(frame, design_spec(~ 0 + mv(prefix)))

  expect_identical(term_data(compiled)$component_id, c("x", "x1"))
  expect_identical(term_data(compiled)$source_block, c("prefix", "prefix"))
  expect_identical(term_data(compiled)$alignment, c("observation", "observation"))
})

test_that("mv syntax and component selectors are strict and deterministic", {
  frame <- .design_frame_fixture()

  expect_error(compile_design(frame, design_spec(~ mv())), "block")
  expect_error(compile_design(frame, design_spec(~ mv(motion, 1, 2))), "arguments")
  expect_error(compile_design(frame, design_spec(~ mv("motion"))), "unquoted block")
  expect_error(compile_design(frame, design_spec(~ mv(motion + age))), "unquoted block")
  expect_error(compile_design(frame, design_spec(~ mv(motion, c(1, 1)))), "unique")
  expect_error(compile_design(frame, design_spec(~ mv(motion, 1.5))), "integers")
  expect_error(compile_design(frame, design_spec(~ mv(motion, TRUE))), "IDs or integer")
  expect_error(
    compile_design(frame, design_spec(~ scale(mv(motion)))),
    "formula algebra"
  )

  named <- compile_design(
    frame,
    design_spec(~ 0 + mv(components = "rotation", block = motion))
  )
  expect_identical(term_data(named)$component_id, "rotation")
})

test_that("mv terms in random formulas fail explicitly until random blocks compile", {
  frame <- .design_frame_fixture()
  expect_error(
    compile_design(
      frame,
      design_spec(~ Fac1, random = ~ 1 + mv(motion) | subject_id)
    ),
    "random-effects"
  )
})

test_that("component metadata remains exact through interactions", {
  frame <- .design_frame_fixture()
  frame$observations$blocks$prefix <- fmridataset::axis_block(
    matrix(seq_len(16L), nrow = 8L),
    components = data.frame(.component_id = c("x", "x1"))
  )
  compiled <- compile_design(frame, design_spec(~ mv(prefix) * Fac1))
  terms <- term_data(compiled)
  interaction <- grepl(":", terms$term, fixed = TRUE)

  expect_identical(
    terms$component_id[terms$term == "prefix[x]"],
    "x"
  )
  expect_identical(
    terms$component_id[terms$term == "prefix[x1]"],
    "x1"
  )
  expect_setequal(
    terms$component_id[interaction],
    c("x", "x1")
  )
  expect_true(all(terms$alignment[interaction] == "observation"))
})

test_that("qualified mv calls compile and tensor blocks fail explicitly", {
  frame <- .design_frame_fixture()
  qualified <- compile_design(
    frame,
    design_spec(~ 0 + multidesign::mv(motion, "rotation"))
  )
  expect_identical(term_data(qualified)$component_id, "rotation")

  tensor <- array(seq_len(8L * 2L * 3L), dim = c(8L, 2L, 3L))
  frame$observations$blocks$tensor <- fmridataset::axis_block(
    tensor,
    components = data.frame(.component_id = c("a", "b"))
  )
  expect_error(
    compile_design(frame, design_spec(~ mv(tensor))),
    "two-dimensional"
  )
})

test_that("random formulas validate every referenced scalar", {
  frame <- .design_frame_fixture()
  expect_error(
    compile_design(
      frame,
      design_spec(~ Fac1, random = ~ 1 + absent_slope | subject_id)
    ),
    "absent_slope"
  )
  expect_error(
    compile_design(frame, design_spec(~ Fac1, random = ~ subject_id)),
    "grouping term"
  )
})

test_that("design compilation is equivariant to row order and feature selection", {
  frame <- .design_frame_fixture()
  spec <- design_spec(~ Fac1 + subject.age + mv(stimulus.visual_pca, 1:2))
  reference <- compile_design(frame, spec)
  order <- c(8L, 2L, 6L, 1L, 7L, 3L, 5L, 4L)
  reordered <- compile_design(frame[order, c(3L, 1L)], spec)

  expect_equal(
    unname(model_matrix(reordered)),
    unname(model_matrix(reference)[order, , drop = FALSE]),
    tolerance = 0,
    ignore_attr = TRUE
  )
  expect_identical(
    reordered$observation_ids,
    reference$observation_ids[order]
  )
  expect_identical(term_data(reordered), term_data(reference))
})
