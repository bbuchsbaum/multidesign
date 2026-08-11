.formula_signature <- function(x) {
  if (is.null(x)) return(NULL)
  paste(deparse(x, width.cutoff = 500L), collapse = "")
}

.design_spec_signature <- function(spec) {
  list(
    schema_version = spec$schema_version,
    fixed = .formula_signature(spec$fixed),
    random = .formula_signature(spec$random),
    contrasts = spec$contrasts,
    na_action = spec$na_action
  )
}

.block_digest_descriptor <- function(block) {
  value <- fmridataset::axis_block_data(block)
  list(
    fingerprint = if (inherits(value, "array_source")) {
      fmridataset::source_fingerprint(value)
    } else {
      digest::digest(value, algo = "sha256", serialize = TRUE)
    },
    component_ids = fmridataset::block_components(block)$.component_id,
    role = block$role,
    units = block$units,
    alignment = if (is.null(block$metadata$.fmridataset_lift$entity)) {
      "observation"
    } else {
      block$metadata$.fmridataset_lift$entity
    }
  )
}

#' Digest semantic inputs to frame-native design compilation
#'
#' The digest includes formulas and policies, stable observation IDs, resolved
#' scalar annotations, and referenced multivariate block fingerprints and
#' component IDs. Imaging assays and the feature axis are intentionally absent.
#'
#' @param frame An `fmri_frame` or synchronized `fmri_view`.
#' @param spec A `design_spec`.
#' @return A lowercase SHA-256 digest.
#' @export
design_input_digest <- function(frame, spec) {
  if (!inherits(frame, c("fmri_frame", "fmri_view"))) {
    stop("`frame` must be an fmri_frame or fmri_view.", call. = FALSE)
  }
  .validate_design_spec(spec)
  .validate_mv_formula_context(spec$fixed[[2L]])
  calls <- .collect_mv_calls(spec$fixed[[2L]])
  keys <- vapply(calls, .mv_call_key, character(1))
  calls <- calls[!duplicated(keys)]
  keys <- keys[!duplicated(keys)]
  parsed <- lapply(calls, .parse_mv_call)
  targets <- vapply(parsed, `[[`, character(1), "target")
  blocks <- fmridataset::obs_blocks(frame, resolve = TRUE)
  missing_blocks <- setdiff(targets, names(blocks))
  if (length(missing_blocks)) {
    stop("Unknown multivariate block: ", missing_blocks[[1L]], call. = FALSE)
  }
  block_descriptors <- lapply(seq_along(calls), function(index) {
    block <- blocks[[targets[[index]]]]
    descriptor <- .block_digest_descriptor(block)
    selected <- .select_block_components(
      block,
      parsed[[index]]$selection,
      environment(spec$fixed)
    )
    descriptor$selected_component_ids <- descriptor$component_ids[selected]
    descriptor$call <- keys[[index]]
    descriptor
  })
  names(block_descriptors) <- keys
  resolved_observations <- as.data.frame(
    fmridataset::observations(frame, resolve = TRUE)
  )
  scalar_variables <- intersect(
    unique(c(all.vars(spec$fixed), all.vars(spec$random))),
    names(resolved_observations)
  )
  digest::digest(
    list(
      schema_version = 1L,
      spec = .design_spec_signature(spec),
      observation_ids = fmridataset::observation_ids(frame),
      observations = resolved_observations[scalar_variables],
      blocks = block_descriptors
    ),
    algo = "sha256",
    serialize = TRUE
  )
}

#' Digest a compiled design
#'
#' @param x A `compiled_design`.
#' @return A lowercase SHA-256 digest independent of formula environments and
#'   runtime cache state.
#' @export
design_digest <- function(x) {
  if (!inherits(x, "compiled_design")) {
    stop("`x` must be a compiled_design.", call. = FALSE)
  }
  digest::digest(
    list(
      schema_version = 1L,
      spec = .design_spec_signature(x$spec),
      model_matrix = x$model_matrix,
      terms = x$terms,
      components = x$components,
      grouping = x$grouping,
      grouping_terms = x$grouping_terms,
      random_effects = x$random_effects,
      observation_ids = x$observation_ids,
      row_audit = x$row_audit
    ),
    algo = "sha256",
    serialize = TRUE
  )
}

.validate_design_cache <- function(cache) {
  if (!inherits(cache, "design_cache") || !is.environment(cache) ||
      !is.environment(cache$values) || !is.numeric(cache$max_entries) ||
      length(cache$max_entries) != 1L) {
    stop("`cache` must be a design_cache.", call. = FALSE)
  }
  invisible(cache)
}

#' Create and inspect an explicit compiled-design cache
#'
#' Cached values are serialized before storage and unserialized on every hit,
#' preventing caller mutation from poisoning later results. The cache is a
#' runtime object and is never embedded in a frame or compiled design.
#'
#' @param max_entries Positive maximum number of compiled designs retained.
#' @param cache A `design_cache`.
#' @return `design_cache()` returns a cache environment; `design_cache_info()`
#'   returns scalar cache statistics; `clear_design_cache()` returns `cache`
#'   invisibly.
#' @export
design_cache <- function(max_entries = 64L) {
  if (!is.numeric(max_entries) || length(max_entries) != 1L ||
      is.na(max_entries) || max_entries != as.integer(max_entries) ||
      max_entries < 1L) {
    stop("`max_entries` must be one positive integer.", call. = FALSE)
  }
  cache <- new.env(parent = emptyenv())
  cache$values <- new.env(hash = TRUE, parent = emptyenv())
  cache$order <- character()
  cache$max_entries <- as.integer(max_entries)
  cache$hits <- 0L
  cache$misses <- 0L
  class(cache) <- "design_cache"
  cache
}

.design_cache_get <- function(cache, key) {
  .validate_design_cache(cache)
  if (!exists(key, envir = cache$values, inherits = FALSE)) {
    cache$misses <- cache$misses + 1L
    return(NULL)
  }
  cache$hits <- cache$hits + 1L
  cache$order <- c(setdiff(cache$order, key), key)
  unserialize(get(key, envir = cache$values, inherits = FALSE))
}

.design_cache_put <- function(cache, key, value) {
  .validate_design_cache(cache)
  assign(key, serialize(value, NULL), envir = cache$values)
  cache$order <- c(setdiff(cache$order, key), key)
  while (length(cache$order) > cache$max_entries) {
    remove(list = cache$order[[1L]], envir = cache$values)
    cache$order <- cache$order[-1L]
  }
  invisible(value)
}

#' @rdname design_cache
#' @export
design_cache_info <- function(cache) {
  .validate_design_cache(cache)
  list(
    entries = as.integer(length(cache$order)),
    hits = cache$hits,
    misses = cache$misses,
    max_entries = cache$max_entries,
    keys = cache$order
  )
}

#' @rdname design_cache
#' @export
clear_design_cache <- function(cache) {
  .validate_design_cache(cache)
  remove(list = ls(cache$values, all.names = TRUE), envir = cache$values)
  cache$order <- character()
  cache$hits <- 0L
  cache$misses <- 0L
  invisible(cache)
}

.normalize_assessment_index <- function(value, ids, label) {
  if (is.character(value)) {
    unknown <- setdiff(value, ids)
    if (length(unknown)) {
      stop(label, " contains unknown observation IDs: ",
        paste(unknown, collapse = ", "), call. = FALSE
      )
    }
    value <- match(value, ids)
  } else if (is.logical(value)) {
    if (length(value) != length(ids) || anyNA(value)) {
      stop(label, " logical selector must be complete and match frame rows.",
        call. = FALSE
      )
    }
    value <- which(value)
  } else if (is.numeric(value)) {
    if (anyNA(value) || any(!is.finite(value)) ||
        any(value != as.integer(value))) {
      stop(label, " must contain finite integer positions.", call. = FALSE)
    }
    value <- as.integer(value)
  } else {
    stop(label, " must contain positions, IDs, or a logical selector.",
      call. = FALSE
    )
  }
  if (!length(value)) stop(label, " must be non-empty.", call. = FALSE)
  if (anyDuplicated(value)) stop(label, " positions must be unique.", call. = FALSE)
  if (any(value < 1L | value > length(ids))) {
    stop(label, " contains positions outside the frame.", call. = FALSE)
  }
  if (length(value) == length(ids)) {
    stop(label, " must leave at least one analysis observation.", call. = FALSE)
  }
  value
}

#' Compile leakage-safe designs for explicit observation folds
#'
#' Each analysis design is compiled independently. Its frozen blueprint is then
#' applied to the assessment view, so factor coding and transformations are
#' learned only from analysis observations.
#'
#' @param frame An `fmri_frame` or synchronized `fmri_view`.
#' @param spec A `design_spec`.
#' @param assessment Non-empty list of assessment selectors. Each selector may
#'   contain integer positions, stable observation IDs, or one logical vector.
#' @param cache Optional `design_cache` used for analysis compilations.
#' @return A `compiled_design_folds` list.
#' @export
compile_design_folds <- function(frame, spec, assessment, cache = NULL) {
  if (!inherits(frame, c("fmri_frame", "fmri_view"))) {
    stop("`frame` must be an fmri_frame or fmri_view.", call. = FALSE)
  }
  .validate_design_spec(spec)
  if (!is.list(assessment) || !length(assessment)) {
    stop("`assessment` must be a non-empty list of fold selectors.", call. = FALSE)
  }
  if (!is.null(cache)) .validate_design_cache(cache)
  ids <- fmridataset::observation_ids(frame)
  indices <- lapply(seq_along(assessment), function(index) {
    .normalize_assessment_index(
      assessment[[index]], ids, paste0("assessment[[", index, "]]"))
  })
  folds <- lapply(seq_along(indices), function(index) {
    assessment_index <- indices[[index]]
    analysis_index <- setdiff(seq_along(ids), assessment_index)
    analysis <- frame[analysis_index, , drop = FALSE]
    assessment_frame <- frame[assessment_index, , drop = FALSE]
    analysis_design <- compile_design(analysis, spec, cache = cache)
    assessment_design <- apply_design(analysis_design, assessment_frame)
    list(
      analysis = analysis,
      assessment = assessment_frame,
      analysis_design = analysis_design,
      assessment_design = assessment_design
    )
  })
  fold_names <- names(assessment)
  if (is.null(fold_names)) fold_names <- paste0("fold_", seq_along(folds))
  names(folds) <- fold_names
  structure(folds, class = c("compiled_design_folds", "list"))
}

#' Observation membership for compiled design folds
#'
#' @param x A `compiled_design_folds` object.
#' @return A normalized data frame with one row per fold, observation, and role.
#' @export
design_fold_data <- function(x) {
  if (!inherits(x, "compiled_design_folds")) {
    stop("`x` must be compiled_design_folds.", call. = FALSE)
  }
  do.call(rbind, lapply(seq_along(x), function(index) {
    data.frame(
      .fold = rep(index, nrow(x[[index]]$analysis) + nrow(x[[index]]$assessment)),
      .obs_id = c(
        fmridataset::observation_ids(x[[index]]$analysis),
        fmridataset::observation_ids(x[[index]]$assessment)
      ),
      role = c(
        rep("analysis", nrow(x[[index]]$analysis)),
        rep("assessment", nrow(x[[index]]$assessment))
      ),
      stringsAsFactors = FALSE
    )
  }))
}

#' @export
print.compiled_design_folds <- function(x, ...) {
  cat("<compiled_design_folds>", length(x), "folds\n")
  invisible(x)
}
