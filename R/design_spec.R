#' Define a frame-native model design
#'
#' A design specification preserves formulas and their environments. It does
#' not freeze a model matrix: coding, missing-value handling, and multivariate
#' block expansion occur when the specification is compiled against a frame.
#'
#' @param fixed One-sided fixed-effects formula.
#' @param random Optional one-sided random-effects formula.
#' @param contrasts Optional named contrast list passed to
#'   [stats::model.matrix()].
#' @param na_action Missing-value policy, either `"fail"` or `"omit"`.
#' @return A serializable `design_spec`.
#' @export
design_spec <- function(
  fixed,
  random = NULL,
  contrasts = NULL,
  na_action = c("fail", "omit")
) {
  if (!inherits(fixed, "formula") || length(fixed) != 2L) {
    stop("`fixed` must be a one-sided formula.", call. = FALSE)
  }
  if (!is.null(random) &&
    (!inherits(random, "formula") || length(random) != 2L)) {
    stop("`random` must be NULL or a one-sided formula.", call. = FALSE)
  }
  if (!is.null(contrasts) &&
    (!is.list(contrasts) || is.null(names(contrasts)))) {
    stop("`contrasts` must be NULL or a named list.", call. = FALSE)
  }
  structure(
    list(
      fixed = fixed,
      random = random,
      contrasts = contrasts,
      na_action = match.arg(na_action),
      schema_version = 1L
    ),
    class = "design_spec"
  )
}

#' Multivariate block formula special
#'
#' `mv()` marks a named observation- or entity-aligned `axis_block` for
#' expansion by [compile_design()]. It is meaningful only inside a formula.
#'
#' @param block Unquoted block name, optionally qualified by an entity name,
#'   such as `motion` or `stimulus.visual_pca`.
#' @param components Optional component positions or stable component IDs.
#' @return This function always errors when evaluated directly.
#' @export
mv <- function(block, components = NULL) {
  stop("`mv()` is a formula special and must be used inside `design_spec()`.",
    call. = FALSE
  )
}

.is_mv_call <- function(x) {
  is.call(x) && identical(x[[1L]], as.name("mv"))
}

.mv_call_key <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")

.collect_mv_calls <- function(x) {
  if (.is_mv_call(x)) {
    return(list(x))
  }
  if (!is.call(x)) {
    return(list())
  }
  unlist(lapply(as.list(x)[-1L], .collect_mv_calls), recursive = FALSE)
}

.replace_mv_calls <- function(x, replacements) {
  if (.is_mv_call(x)) {
    return(replacements[[.mv_call_key(x)]])
  }
  if (!is.call(x)) {
    return(x)
  }
  as.call(lapply(as.list(x), .replace_mv_calls, replacements = replacements))
}

.plus_expression <- function(names) {
  expressions <- lapply(names, as.name)
  if (length(expressions) == 1L) {
    return(expressions[[1L]])
  }
  Reduce(function(left, right) call("+", left, right), expressions)
}

.frame_entities <- function(frame) {
  if (inherits(frame, "fmri_view")) frame$base$entities else frame$entities
}

.block_matrix <- function(block) {
  value <- fmridataset::axis_block_data(block)
  if (inherits(value, "array_source")) {
    value <- fmridataset::source_read(value)
  }
  as.matrix(value)
}

.select_block_components <- function(block, selection, environment) {
  components <- fmridataset::block_components(block)
  ids <- components$.component_id
  if (is.null(selection)) {
    return(seq_along(ids))
  }
  selection <- eval(selection, envir = environment)
  if (is.character(selection)) {
    index <- match(selection, ids)
  } else {
    index <- as.integer(selection)
  }
  if (!length(index) || anyNA(index) || any(index < 1L | index > length(ids))) {
    stop("`mv()` component selection is invalid for this block.", call. = FALSE)
  }
  index
}

.entity_block <- function(frame, entity_name, block_name, observations) {
  entity <- .frame_entities(frame)[[entity_name]]
  if (is.null(entity) || is.null(entity$blocks[[block_name]])) {
    stop("Unknown entity block: ", entity_name, ".", block_name, call. = FALSE)
  }
  entity_data <- entity$data
  if (!is.data.frame(entity_data)) {
    stop("Entity block requires an entity data frame.", call. = FALSE)
  }
  key <- entity$key
  if (is.null(key)) {
    preferred <- paste0(entity_name, "_id")
    candidates <- intersect(names(observations), names(entity_data))
    key <- if (preferred %in% candidates) preferred else candidates[1L]
  }
  if (length(key) != 1L || is.na(key) || !nzchar(key)) {
    stop("Could not resolve the observation-to-entity key for ", entity_name, ".",
      call. = FALSE
    )
  }
  index <- match(observations[[key]], entity_data[[key]])
  if (anyNA(index)) {
    stop("Entity relation contains unresolved observation keys.", call. = FALSE)
  }
  list(block = entity$blocks[[block_name]], index = index)
}

.resolve_mv_block <- function(frame, call, observations, environment) {
  target <- paste(deparse(call[[2L]], width.cutoff = 500L), collapse = "")
  observation_blocks <- fmridataset::obs_blocks(frame)
  alignment <- "observation"
  row_index <- seq_len(nrow(observations))

  if (target %in% names(observation_blocks)) {
    block <- observation_blocks[[target]]
  } else {
    entities <- names(.frame_entities(frame))
    matching <- entities[startsWith(target, paste0(entities, "."))]
    if (!length(matching)) {
      stop("Unknown multivariate block: ", target, call. = FALSE)
    }
    entity_name <- matching[[which.max(nchar(matching))]]
    block_name <- substring(target, nchar(entity_name) + 2L)
    resolved <- .entity_block(frame, entity_name, block_name, observations)
    block <- resolved$block
    row_index <- resolved$index
    alignment <- entity_name
  }

  selection <- if (length(call) >= 3L) call[[3L]] else NULL
  component_index <- .select_block_components(
    block,
    selection,
    environment
  )
  components <- fmridataset::block_components(block)[component_index, , drop = FALSE]
  matrix <- .block_matrix(block)[row_index, component_index, drop = FALSE]
  list(
    target = target,
    matrix = matrix,
    components = components,
    alignment = alignment
  )
}

.random_group_vars <- function(formula) {
  if (is.null(formula)) {
    return(character())
  }
  find_groups <- function(x) {
    if (!is.call(x)) {
      return(character())
    }
    if (identical(x[[1L]], as.name("|")) || identical(x[[1L]], as.name("||"))) {
      return(all.vars(x[[3L]]))
    }
    unique(unlist(lapply(as.list(x)[-1L], find_groups), use.names = FALSE))
  }
  find_groups(formula[[2L]])
}

.human_term <- function(term, component_map) {
  out <- term
  for (index in seq_len(nrow(component_map))) {
    label <- paste0(
      component_map$block[[index]],
      "[",
      component_map$component_id[[index]],
      "]"
    )
    out <- gsub(component_map$column[[index]], label, out, fixed = TRUE)
  }
  out
}

#' Compile a design specification against an fmri frame
#'
#' @param frame An `fmri_frame` or synchronized `fmri_view`.
#' @param spec A `design_spec`.
#' @return A `compiled_design` containing the dense fixed-effects model matrix,
#'   component-aware term metadata, grouping data, and retained source spec.
#' @export
compile_design <- function(frame, spec) {
  if (!inherits(frame, c("fmri_frame", "fmri_view"))) {
    stop("`frame` must be an fmri_frame or fmri_view.", call. = FALSE)
  }
  if (!inherits(spec, "design_spec")) {
    stop("`spec` must be a design_spec.", call. = FALSE)
  }
  data <- as.data.frame(fmridataset::observations(frame))
  calls <- .collect_mv_calls(spec$fixed[[2L]])
  keys <- vapply(calls, .mv_call_key, character(1))
  calls <- calls[!duplicated(keys)]
  keys <- keys[!duplicated(keys)]
  replacements <- list()
  component_rows <- list()

  for (index in seq_along(calls)) {
    resolved <- .resolve_mv_block(
      frame,
      calls[[index]],
      data,
      environment(spec$fixed)
    )
    component_ids <- resolved$components$.component_id
    columns <- make.names(
      paste("mv", resolved$target, component_ids, sep = "__"),
      unique = TRUE
    )
    if (any(columns %in% names(data))) {
      stop("Compiled multivariate column name collides with scalar metadata.",
        call. = FALSE
      )
    }
    data[columns] <- resolved$matrix
    replacements[[keys[[index]]]] <- .plus_expression(columns)
    component_rows[[index]] <- data.frame(
      column = columns,
      block = resolved$target,
      component_id = component_ids,
      alignment = resolved$alignment,
      stringsAsFactors = FALSE
    )
  }
  component_map <- if (length(component_rows)) {
    do.call(rbind, component_rows)
  } else {
    data.frame(
      column = character(),
      block = character(),
      component_id = character(),
      alignment = character()
    )
  }

  rewritten <- spec$fixed
  rewritten[[2L]] <- .replace_mv_calls(rewritten[[2L]], replacements)
  environment(rewritten) <- environment(spec$fixed)
  na_function <- if (identical(spec$na_action, "fail")) stats::na.fail else stats::na.omit
  model_frame <- stats::model.frame(rewritten, data = data, na.action = na_function)
  matrix <- stats::model.matrix(
    rewritten,
    data = model_frame,
    contrasts.arg = spec$contrasts
  )
  terms <- stats::terms(rewritten, data = model_frame)
  term_labels <- c("(Intercept)", attr(terms, "term.labels"))
  assigned <- attr(matrix, "assign") + 1L

  term_table <- lapply(seq_len(ncol(matrix)), function(column_index) {
    model_column <- colnames(matrix)[[column_index]]
    matches <- which(vapply(
      component_map$column,
      function(token) grepl(token, model_column, fixed = TRUE),
      logical(1)
    ))
    data.frame(
      model_column = model_column,
      term = .human_term(term_labels[[assigned[[column_index]]]], component_map),
      source_type = if (length(matches)) "axis_block" else "scalar",
      source_block = if (length(matches)) {
        paste(unique(component_map$block[matches]), collapse = ";")
      } else {
        NA_character_
      },
      component_id = if (length(matches)) {
        paste(unique(component_map$component_id[matches]), collapse = ";")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  term_table <- do.call(rbind, term_table)

  group_vars <- .random_group_vars(spec$random)
  missing_groups <- setdiff(group_vars, names(data))
  if (length(missing_groups)) {
    stop("Random grouping variables are missing: ",
      paste(missing_groups, collapse = ", "),
      call. = FALSE
    )
  }
  grouping <- data[group_vars]
  retained_rows <- match(rownames(model_frame), rownames(data))
  if (anyNA(retained_rows)) {
    stop("Could not retain observation identity after missing-value handling.",
      call. = FALSE
    )
  }

  structure(
    list(
      model_matrix = matrix,
      terms = term_table,
      grouping = grouping[retained_rows, , drop = FALSE],
      observation_ids = fmridataset::observation_ids(frame)[retained_rows],
      fixed_formula = rewritten,
      random_formula = spec$random,
      component_map = component_map,
      spec = spec
    ),
    class = "compiled_design"
  )
}

#' Access compiled design products
#'
#' @param x A `compiled_design`.
#' @return The requested matrix or metadata table.
#' @export
model_matrix <- function(x) {
  if (!inherits(x, "compiled_design")) stop("`x` must be a compiled_design.", call. = FALSE)
  x$model_matrix
}

#' @rdname model_matrix
#' @export
term_data <- function(x) {
  if (!inherits(x, "compiled_design")) stop("`x` must be a compiled_design.", call. = FALSE)
  x$terms
}

#' @rdname model_matrix
#' @export
grouping_data <- function(x) {
  if (!inherits(x, "compiled_design")) stop("`x` must be a compiled_design.", call. = FALSE)
  x$grouping
}

#' @export
print.compiled_design <- function(x, ...) {
  cat(
    "<compiled_design>", nrow(x$model_matrix), "observations x",
    ncol(x$model_matrix), "columns\n"
  )
  cat("  multivariate components:", nrow(x$component_map), "\n")
  cat(
    "  grouping variables:",
    if (ncol(x$grouping)) paste(names(x$grouping), collapse = ", ") else "none",
    "\n"
  )
  invisible(x)
}
