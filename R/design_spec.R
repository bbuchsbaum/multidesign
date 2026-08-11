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
#' @details Fixed formulas are compiled against relation-resolved observation
#'   metadata. Random formulas must contain at least one grouping term; `mv()`
#'   is currently restricted to fixed-effects formula algebra.
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
  if (!is.null(contrasts) && !is.list(contrasts)) {
    stop("`contrasts` must be NULL or a named list.", call. = FALSE)
  }
  if (!is.null(contrasts)) {
    contrast_names <- names(contrasts)
    if (is.null(contrast_names) || anyNA(contrast_names) ||
        any(!nzchar(contrast_names))) {
      stop("`contrasts` names must be non-empty strings.", call. = FALSE)
    }
    if (anyDuplicated(contrast_names)) {
      stop("`contrasts` names must be unique.", call. = FALSE)
    }
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
#' expansion by [compile_design()]. Blocks are resolved through the frame's
#' relation registry and must be two dimensional. It is meaningful only as a
#' term in fixed-effects formula algebra.
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
  if (!is.call(x)) return(FALSE)
  if (identical(x[[1L]], as.name("mv"))) return(TRUE)
  head <- x[[1L]]
  is.call(head) && identical(head[[1L]], as.name("::")) &&
    identical(as.character(head[[2L]]), "multidesign") &&
    identical(as.character(head[[3L]]), "mv")
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

.validate_mv_formula_context <- function(x) {
  if (.is_mv_call(x) || !is.call(x)) return(invisible(TRUE))
  contains_mv <- length(.collect_mv_calls(x)) > 0L
  head <- x[[1L]]
  allowed <- is.symbol(head) && as.character(head) %in% c(
    "+", "-", "*", ":", "/", "^", "("
  )
  if (contains_mv && !allowed) {
    stop(
      "`mv()` must appear as a formula term or within formula algebra; wrapping it in another function is ambiguous.",
      call. = FALSE
    )
  }
  for (value in as.list(x)[-1L]) .validate_mv_formula_context(value)
  invisible(TRUE)
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

.block_matrix <- function(block, component_index) {
  value <- fmridataset::axis_block_data(block)
  if (inherits(value, "array_source")) {
    return(fmridataset::source_read(
      value,
      observations = seq_len(fmridataset::source_shape(value)[[1L]]),
      features = component_index
    ))
  }
  as.matrix(value)[, component_index, drop = FALSE]
}

.parse_mv_call <- function(x) {
  matched <- tryCatch(
    match.call(definition = mv, call = x, expand.dots = FALSE),
    error = function(error) {
      stop(
        "`mv()` accepts only `block` and optional `components` arguments.",
        call. = FALSE
      )
    }
  )
  if (is.null(matched$block)) {
    stop("`mv()` requires an unquoted block name.", call. = FALSE)
  }
  if (!is.symbol(matched$block)) {
    stop(
      "`mv()` block must be one unquoted block name such as `motion` or `stimulus.visual_pca`.",
      call. = FALSE
    )
  }
  list(
    target = as.character(matched$block),
    selection = if (!is.null(matched$components)) matched$components else NULL
  )
}

.select_block_components <- function(block, selection, environment) {
  components <- fmridataset::block_components(block)
  ids <- components$.component_id
  if (is.null(selection)) {
    return(seq_along(ids))
  }
  selection <- eval(selection, envir = environment)
  if (is.null(selection)) return(seq_along(ids))
  if (is.character(selection)) {
    if (!length(selection) || anyNA(selection) || any(!nzchar(selection)) ||
        anyDuplicated(selection)) {
      stop("`mv()` component IDs must be unique non-empty strings.", call. = FALSE)
    }
    index <- match(selection, ids)
  } else if (is.numeric(selection) && !is.logical(selection)) {
    if (!length(selection) || anyNA(selection) || any(!is.finite(selection)) ||
        any(selection != as.integer(selection))) {
      stop("`mv()` component positions must be finite integers.", call. = FALSE)
    }
    index <- as.integer(selection)
    if (anyDuplicated(index)) {
      stop("`mv()` component positions must be unique.", call. = FALSE)
    }
  } else {
    stop("`mv()` components must be stable component IDs or integer positions.",
      call. = FALSE
    )
  }
  if (!length(index) || anyNA(index) || any(index < 1L | index > length(ids))) {
    stop("`mv()` component selection is invalid for this block.", call. = FALSE)
  }
  index
}

.resolve_mv_block <- function(frame, call, observations, environment) {
  parsed <- .parse_mv_call(call)
  target <- parsed$target
  observation_blocks <- fmridataset::obs_blocks(frame, resolve = TRUE)
  alignment <- "observation"
  if (!target %in% names(observation_blocks)) {
    stop("Unknown multivariate block: ", target, call. = FALSE)
  }
  block <- observation_blocks[[target]]
  block_data <- fmridataset::axis_block_data(block)
  block_shape <- if (inherits(block_data, "array_source")) {
    fmridataset::source_shape(block_data)
  } else {
    dim(block_data)
  }
  if (length(block_shape) != 2L) {
    stop("`mv()` requires a two-dimensional axis block.", call. = FALSE)
  }
  lift <- block$metadata$.fmridataset_lift
  if (is.list(lift) && is.character(lift$entity) && length(lift$entity) == 1L) {
    alignment <- lift$entity
  }

  component_index <- .select_block_components(
    block,
    parsed$selection,
    environment
  )
  components <- fmridataset::block_components(block)[component_index, , drop = FALSE]
  matrix <- .block_matrix(block, component_index)
  if (nrow(matrix) != nrow(observations)) {
    stop("Resolved multivariate block is not aligned to observations.", call. = FALSE)
  }
  list(
    target = target,
    matrix = matrix,
    components = components,
    alignment = alignment
  )
}

.resolve_frozen_mv_block <- function(frame, term, observations) {
  observation_blocks <- fmridataset::obs_blocks(frame, resolve = TRUE)
  if (!term$target %in% names(observation_blocks)) {
    stop("Unknown multivariate block: ", term$target, call. = FALSE)
  }
  block <- observation_blocks[[term$target]]
  components <- fmridataset::block_components(block)
  component_index <- match(term$component_ids, components$.component_id)
  if (anyNA(component_index)) {
    stop(
      "Multivariate block `", term$target,
      "` is missing frozen component IDs: ",
      paste(term$component_ids[is.na(component_index)], collapse = ", "),
      call. = FALSE
    )
  }
  matrix <- .block_matrix(block, component_index)
  if (length(dim(matrix)) != 2L || nrow(matrix) != nrow(observations)) {
    stop("Resolved multivariate block is not aligned to observations.", call. = FALSE)
  }
  list(matrix = matrix, components = components[component_index, , drop = FALSE])
}

.random_group_vars <- function(formula) {
  if (is.null(formula)) {
    return(character())
  }
  calls <- .collect_random_calls(formula[[2L]])
  unique(unlist(lapply(calls, function(x) all.vars(x[[3L]])), use.names = FALSE))
}

.collect_random_calls <- function(x) {
  if (!is.call(x)) return(list())
  if (identical(x[[1L]], as.name("|")) || identical(x[[1L]], as.name("||"))) {
    return(list(x))
  }
  unlist(lapply(as.list(x)[-1L], .collect_random_calls), recursive = FALSE)
}

.human_term <- function(term, component_map) {
  out <- term
  order <- order(nchar(component_map$column), decreasing = TRUE)
  for (index in order) {
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

.term_component_indices <- function(term, component_map) {
  if (!nrow(component_map) || identical(term, "(Intercept)")) return(integer())
  variables <- tryCatch(all.vars(str2lang(term)), error = function(error) character())
  which(component_map$column %in% variables)
}

.validate_design_spec <- function(spec) {
  required <- c(
    "fixed", "random", "contrasts", "na_action", "schema_version"
  )
  if (!inherits(spec, "design_spec") ||
      !identical(names(unclass(spec)), required) ||
      !identical(spec$schema_version, 1L)) {
    stop("`spec` must be a valid design_spec.", call. = FALSE)
  }
  design_spec(
    fixed = spec$fixed,
    random = spec$random,
    contrasts = spec$contrasts,
    na_action = spec$na_action
  )
  invisible(spec)
}

.prepare_design_data <- function(frame, spec, frozen_mv_terms = NULL) {
  data <- as.data.frame(fmridataset::observations(frame, resolve = TRUE))
  calls <- .collect_mv_calls(spec$fixed[[2L]])
  keys <- vapply(calls, .mv_call_key, character(1))
  calls <- calls[!duplicated(keys)]
  keys <- keys[!duplicated(keys)]
  replacements <- list()
  component_rows <- list()
  mv_terms <- vector("list", length(calls))

  if (!is.null(frozen_mv_terms) && length(frozen_mv_terms) != length(calls)) {
    stop("The frozen multivariate design does not match the source formula.",
      call. = FALSE
    )
  }

  for (index in seq_along(calls)) {
    if (is.null(frozen_mv_terms)) {
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
      term <- list(
        key = keys[[index]],
        target = resolved$target,
        component_ids = component_ids,
        columns = columns,
        alignment = resolved$alignment
      )
    } else {
      term <- frozen_mv_terms[[index]]
      if (!identical(term$key, keys[[index]])) {
        stop("The frozen multivariate design does not match the source formula.",
          call. = FALSE
        )
      }
      resolved_block <- .resolve_frozen_mv_block(frame, term, data)
      resolved <- list(
        matrix = resolved_block$matrix,
        components = resolved_block$components,
        target = term$target,
        alignment = term$alignment
      )
      component_ids <- term$component_ids
      columns <- term$columns
    }
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
    mv_terms[[index]] <- term
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
  list(
    data = data,
    rewritten = rewritten,
    component_map = component_map,
    mv_terms = mv_terms
  )
}

.design_missingness <- function(data, fixed_vars, random_vars) {
  missing_for <- function(variables) {
    if (!length(variables)) return(rep(FALSE, nrow(data)))
    !stats::complete.cases(data[, variables, drop = FALSE])
  }
  list(
    fixed = missing_for(fixed_vars),
    random = missing_for(random_vars)
  )
}

.random_complete_cases <- function(data, formula) {
  valid <- rep(TRUE, nrow(data))
  if (is.null(formula) || !nrow(data)) return(valid)
  formula_environment <- environment(formula)
  for (random_call in .collect_random_calls(formula[[2L]])) {
    slope_formula <- call("~", random_call[[2L]])
    environment(slope_formula) <- formula_environment
    model_frame <- stats::model.frame(
      slope_formula,
      data = data,
      na.action = stats::na.pass
    )
    grouping_value <- eval(
      random_call[[3L]],
      envir = data,
      enclos = formula_environment
    )
    if (is.data.frame(grouping_value) || is.matrix(grouping_value) ||
        length(grouping_value) != nrow(data)) {
      stop("Random-effects grouping expressions must produce one value per observation.",
        call. = FALSE
      )
    }
    valid <- valid & stats::complete.cases(model_frame) & !is.na(grouping_value)
  }
  valid
}

.missing_design_error <- function(ids, missing) {
  affected <- ids[missing]
  stop(
    length(affected), " observation", if (length(affected) == 1L) "" else "s",
    " contain missing values in design inputs: ", paste(affected, collapse = ", "),
    call. = FALSE
  )
}

.validate_contrast_predictors <- function(contrasts, model_frame) {
  if (is.null(contrasts)) return(invisible(TRUE))
  invalid <- names(contrasts)[
    !names(contrasts) %in% names(model_frame) |
      !vapply(names(contrasts), function(name) {
        name %in% names(model_frame) && is.factor(model_frame[[name]])
      }, logical(1))
  ]
  if (length(invalid)) {
    stop(
      "Contrast declaration is unused or not a categorical predictor: ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.freeze_contrasts <- function(model_frame, matrix, requested) {
  schemes <- attr(matrix, "contrasts")
  if (is.null(schemes)) return(NULL)
  out <- lapply(names(schemes), function(name) {
    value <- model_frame[[name]]
    scheme <- if (name %in% names(requested)) requested[[name]] else schemes[[name]]
    stats::contrasts(value) <- scheme
    frozen <- stats::contrasts(value)
    rownames(frozen) <- levels(value)
    frozen
  })
  names(out) <- names(schemes)
  out
}

.empty_component_data <- function() {
  data.frame(
    model_column = character(),
    column_index = integer(),
    term_id = character(),
    generated_column = character(),
    source_block = character(),
    component_id = character(),
    alignment = character(),
    stringsAsFactors = FALSE
  )
}

.empty_term_data <- function() {
  data.frame(
    model_column = character(),
    column_index = integer(),
    term_id = character(),
    term_index = integer(),
    term = character(),
    source_term = character(),
    is_intercept = logical(),
    source_type = character(),
    component_count = integer(),
    source_block = character(),
    component_id = character(),
    alignment = character(),
    stringsAsFactors = FALSE
  )
}

.build_term_metadata <- function(matrix, terms, component_map) {
  term_labels <- c("(Intercept)", attr(terms, "term.labels"))
  term_indices <- attr(matrix, "assign")
  assigned <- term_indices + 1L
  component_rows <- list()
  rows <- lapply(seq_len(ncol(matrix)), function(column_index) {
    model_column <- colnames(matrix)[[column_index]]
    source_term <- term_labels[[assigned[[column_index]]]]
    matches <- .term_component_indices(source_term, component_map)
    variables <- if (identical(source_term, "(Intercept)")) {
      character()
    } else {
      tryCatch(all.vars(str2lang(source_term)), error = function(error) character())
    }
    scalar_variables <- setdiff(variables, component_map$column[matches])
    source_type <- if (identical(source_term, "(Intercept)")) {
      "intercept"
    } else if (length(matches) && length(scalar_variables)) {
      "mixed"
    } else if (length(matches)) {
      "axis_block"
    } else {
      "scalar"
    }
    term_id <- paste0("fixed:", term_indices[[column_index]])
    if (length(matches)) {
      component_rows[[length(component_rows) + 1L]] <<- data.frame(
        model_column = rep(model_column, length(matches)),
        column_index = rep(column_index, length(matches)),
        term_id = rep(term_id, length(matches)),
        generated_column = component_map$column[matches],
        source_block = component_map$block[matches],
        component_id = component_map$component_id[matches],
        alignment = component_map$alignment[matches],
        stringsAsFactors = FALSE
      )
    }
    data.frame(
      model_column = model_column,
      column_index = column_index,
      term_id = term_id,
      term_index = term_indices[[column_index]],
      term = .human_term(source_term, component_map),
      source_term = source_term,
      is_intercept = identical(source_term, "(Intercept)"),
      source_type = source_type,
      component_count = length(matches),
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
      alignment = if (length(matches)) {
        paste(unique(component_map$alignment[matches]), collapse = ";")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  list(
    terms = if (length(rows)) do.call(rbind, rows) else .empty_term_data(),
    components = if (length(component_rows)) {
      do.call(rbind, component_rows)
    } else {
      .empty_component_data()
    }
  )
}

.empty_random_effect_data <- function() {
  data.frame(
    random_term_id = character(),
    effect_column = character(),
    effect_index = integer(),
    effect_term = character(),
    is_intercept = logical(),
    grouping_expression = character(),
    operator = character(),
    correlated = logical(),
    n_groups = integer(),
    stringsAsFactors = FALSE
  )
}

.empty_grouping_term_data <- function() {
  data.frame(
    random_term_id = character(),
    grouping_expression = character(),
    grouping_variable = character(),
    grouping_variable_index = integer(),
    operator = character(),
    correlated = logical(),
    n_groups = integer(),
    stringsAsFactors = FALSE
  )
}

.build_random_metadata <- function(data, formula, frozen_terms = NULL) {
  if (is.null(formula)) {
    return(list(
      effects = .empty_random_effect_data(),
      groups = .empty_grouping_term_data(),
      blueprints = list()
    ))
  }
  calls <- .collect_random_calls(formula[[2L]])
  if (!is.null(frozen_terms) && length(frozen_terms) != length(calls)) {
    stop("The frozen random-effects design does not match the source formula.",
      call. = FALSE
    )
  }
  effect_rows <- list()
  group_rows <- list()
  blueprints <- vector("list", length(calls))
  formula_environment <- environment(formula)

  for (index in seq_along(calls)) {
    random_call <- calls[[index]]
    operator <- as.character(random_call[[1L]])
    random_term_id <- paste0("random:", index)
    grouping_expression <- paste(deparse(random_call[[3L]], width.cutoff = 500L), collapse = "")
    grouping_variables <- all.vars(random_call[[3L]])
    grouping_value <- eval(random_call[[3L]], envir = data, enclos = formula_environment)
    if (is.data.frame(grouping_value) || is.matrix(grouping_value) ||
        length(grouping_value) != nrow(data)) {
      stop("Random-effects grouping expressions must produce one value per observation.",
        call. = FALSE
      )
    }
    n_groups <- length(unique(grouping_value[!is.na(grouping_value)]))
    slope_formula <- call("~", random_call[[2L]])
    environment(slope_formula) <- formula_environment

    if (is.null(frozen_terms)) {
      model_frame <- stats::model.frame(
        slope_formula,
        data = data,
        na.action = stats::na.fail
      )
      slope_terms <- stats::terms(model_frame)
      matrix <- stats::model.matrix(slope_terms, data = model_frame)
      frozen_contrasts <- .freeze_contrasts(model_frame, matrix, NULL)
      xlevels <- lapply(
        model_frame[vapply(model_frame, is.factor, logical(1))],
        levels
      )
      blueprint <- list(
        random_term_id = random_term_id,
        terms = slope_terms,
        xlevels = xlevels,
        contrasts = frozen_contrasts,
        model_columns = colnames(matrix),
        grouping_expression = grouping_expression,
        grouping_variables = grouping_variables,
        operator = operator
      )
    } else {
      blueprint <- frozen_terms[[index]]
      if (!identical(blueprint$random_term_id, random_term_id) ||
          !identical(blueprint$grouping_expression, grouping_expression) ||
          !identical(blueprint$operator, operator)) {
        stop("The frozen random-effects design does not match the source formula.",
          call. = FALSE
        )
      }
      model_frame <- stats::model.frame(
        blueprint$terms,
        data = data,
        na.action = stats::na.fail,
        xlev = blueprint$xlevels
      )
      slope_terms <- blueprint$terms
      matrix <- stats::model.matrix(
        slope_terms,
        data = model_frame,
        contrasts.arg = blueprint$contrasts
      )
      if (!identical(colnames(matrix), blueprint$model_columns)) {
        stop("Applied random-effect columns do not match the frozen blueprint.",
          call. = FALSE
        )
      }
    }
    blueprints[[index]] <- blueprint
    assigned <- attr(matrix, "assign")
    labels <- c("(Intercept)", attr(slope_terms, "term.labels"))
    effect_rows[[index]] <- if (ncol(matrix)) {
      data.frame(
        random_term_id = rep(random_term_id, ncol(matrix)),
        effect_column = colnames(matrix),
        effect_index = seq_len(ncol(matrix)),
        effect_term = labels[assigned + 1L],
        is_intercept = assigned == 0L,
        grouping_expression = rep(grouping_expression, ncol(matrix)),
        operator = rep(operator, ncol(matrix)),
        correlated = rep(identical(operator, "|"), ncol(matrix)),
        n_groups = rep(as.integer(n_groups), ncol(matrix)),
        stringsAsFactors = FALSE
      )
    } else {
      .empty_random_effect_data()
    }
    group_rows[[index]] <- data.frame(
      random_term_id = rep(random_term_id, length(grouping_variables)),
      grouping_expression = rep(grouping_expression, length(grouping_variables)),
      grouping_variable = grouping_variables,
      grouping_variable_index = seq_along(grouping_variables),
      operator = rep(operator, length(grouping_variables)),
      correlated = rep(identical(operator, "|"), length(grouping_variables)),
      n_groups = rep(as.integer(n_groups), length(grouping_variables)),
      stringsAsFactors = FALSE
    )
  }
  list(
    effects = if (length(effect_rows)) do.call(rbind, effect_rows) else {
      .empty_random_effect_data()
    },
    groups = if (length(group_rows)) do.call(rbind, group_rows) else {
      .empty_grouping_term_data()
    },
    blueprints = blueprints
  )
}

.assemble_compiled_design <- function(
  frame,
  spec,
  prepared,
  frozen_blueprint = NULL
) {
  data <- prepared$data
  group_vars <- .random_group_vars(spec$random)
  random_vars <- if (is.null(spec$random)) character() else all.vars(spec$random)
  fixed_vars <- all.vars(prepared$rewritten)
  missing_vars <- setdiff(unique(c(fixed_vars, random_vars)), names(data))
  if (length(missing_vars)) {
    stop("Design variables are missing: ", paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  missingness <- .design_missingness(data, fixed_vars, random_vars)
  input_missing <- missingness$fixed | missingness$random
  ids <- fmridataset::observation_ids(frame)
  if (identical(spec$na_action, "fail") && any(input_missing)) {
    .missing_design_error(ids, input_missing)
  }
  candidate <- if (identical(spec$na_action, "omit")) !input_missing else {
    rep(TRUE, nrow(data))
  }
  initial_candidates <- which(candidate)
  random_complete <- .random_complete_cases(
    data[initial_candidates, , drop = FALSE],
    spec$random
  )
  random_transform_missing <- rep(FALSE, nrow(data))
  random_transform_missing[initial_candidates[!random_complete]] <- TRUE
  if (identical(spec$na_action, "fail") && any(random_transform_missing)) {
    .missing_design_error(ids, random_transform_missing)
  }
  if (identical(spec$na_action, "omit")) {
    candidate[random_transform_missing] <- FALSE
    missingness$random[random_transform_missing] <- TRUE
  }
  candidate_rows <- which(candidate)
  if (!length(candidate_rows)) {
    stop("No observations remain after missing-value handling.", call. = FALSE)
  }
  candidate_data <- data[candidate_rows, , drop = FALSE]
  na_function <- if (identical(spec$na_action, "fail")) stats::na.fail else stats::na.omit

  if (is.null(frozen_blueprint)) {
    model_frame <- stats::model.frame(
      prepared$rewritten,
      data = candidate_data,
      na.action = na_function
    )
    terms <- stats::terms(model_frame)
    .validate_contrast_predictors(spec$contrasts, model_frame)
    initial_matrix <- stats::model.matrix(
      terms,
      data = model_frame,
      contrasts.arg = spec$contrasts
    )
    frozen_contrasts <- .freeze_contrasts(model_frame, initial_matrix, spec$contrasts)
    matrix <- initial_matrix
    xlevels <- lapply(model_frame[vapply(model_frame, is.factor, logical(1))], levels)
    blueprint <- list(
      schema_version = 1L,
      terms = terms,
      xlevels = xlevels,
      contrasts = frozen_contrasts,
      model_columns = colnames(matrix),
      fixed_vars = fixed_vars,
      random_vars = random_vars,
      group_vars = group_vars,
      mv_terms = prepared$mv_terms
    )
    class(blueprint) <- "design_blueprint"
  } else {
    blueprint <- frozen_blueprint
    model_frame <- stats::model.frame(
      blueprint$terms,
      data = candidate_data,
      na.action = na_function,
      xlev = blueprint$xlevels
    )
    terms <- blueprint$terms
    matrix <- stats::model.matrix(
      terms,
      data = model_frame,
      contrasts.arg = blueprint$contrasts
    )
    if (!identical(colnames(matrix), blueprint$model_columns)) {
      stop("Applied design columns do not match the frozen blueprint.",
        call. = FALSE
      )
    }
  }

  retained_in_candidate <- match(rownames(model_frame), rownames(candidate_data))
  if (anyNA(retained_in_candidate)) {
    stop("Could not retain observation identity after missing-value handling.",
      call. = FALSE
    )
  }
  retained_rows <- candidate_rows[retained_in_candidate]
  retained <- rep(FALSE, nrow(data))
  retained[retained_rows] <- TRUE
  transform_missing <- !retained & candidate
  missingness$fixed[transform_missing] <- TRUE
  row_audit <- data.frame(
    .obs_id = ids,
    fixed_missing = missingness$fixed,
    random_missing = missingness$random,
    retained = retained,
    stringsAsFactors = FALSE
  )
  random_metadata <- .build_random_metadata(
    data[retained_rows, , drop = FALSE],
    spec$random,
    frozen_terms = if (is.null(frozen_blueprint)) {
      NULL
    } else {
      blueprint$random_terms
    }
  )
  if (is.null(frozen_blueprint)) {
    blueprint$random_terms <- random_metadata$blueprints
  }
  fixed_metadata <- .build_term_metadata(
    matrix,
    terms,
    prepared$component_map
  )

  structure(
    list(
      model_matrix = matrix,
      terms = fixed_metadata$terms,
      components = fixed_metadata$components,
      grouping = data[retained_rows, group_vars, drop = FALSE],
      grouping_terms = random_metadata$groups,
      random_effects = random_metadata$effects,
      observation_ids = ids[retained_rows],
      fixed_formula = prepared$rewritten,
      random_formula = spec$random,
      component_map = prepared$component_map,
      row_audit = row_audit,
      blueprint = blueprint,
      spec = spec
    ),
    class = "compiled_design"
  )
}

#' Compile a design specification against an fmri frame
#'
#' @param frame An `fmri_frame` or synchronized `fmri_view`.
#' @param spec A `design_spec`.
#' @param cache Optional runtime [design_cache()]. Cache keys exclude imaging
#'   assays and feature layout and include every semantic design dependency.
#' @return A `compiled_design` containing the dense fixed-effects model matrix,
#'   exact component- and alignment-aware term metadata, grouping data, stable
#'   observation IDs, and the retained source specification.
#' @export
compile_design <- function(frame, spec, cache = NULL) {
  if (!inherits(frame, c("fmri_frame", "fmri_view"))) {
    stop("`frame` must be an fmri_frame or fmri_view.", call. = FALSE)
  }
  if (!inherits(spec, "design_spec")) {
    stop("`spec` must be a design_spec.", call. = FALSE)
  }
  .validate_design_spec(spec)
  if (!is.null(cache)) .validate_design_cache(cache)
  if (!is.null(spec$random) &&
      length(.collect_mv_calls(spec$random[[2L]]))) {
    stop(
      "`mv()` terms in random-effects formulas are not yet supported.",
      call. = FALSE
    )
  }
  .validate_mv_formula_context(spec$fixed[[2L]])
  group_vars <- .random_group_vars(spec$random)
  if (!is.null(spec$random) && !length(group_vars)) {
    stop("`random` must contain at least one `|` or `||` grouping term.",
      call. = FALSE
    )
  }
  input_digest <- design_input_digest(frame, spec)
  if (!is.null(cache)) {
    cached <- .design_cache_get(cache, input_digest)
    if (!is.null(cached)) return(cached)
  }
  prepared <- .prepare_design_data(frame, spec)
  compiled <- .assemble_compiled_design(frame, spec, prepared)
  compiled$input_digest <- input_digest
  if (!is.null(cache)) .design_cache_put(cache, input_digest, compiled)
  compiled
}

#' Apply a frozen design blueprint to another frame
#'
#' Reuses the training factor levels, explicit contrast matrices, formula
#' transformations, model columns, and stable multivariate component IDs.
#'
#' @param x A training `compiled_design`.
#' @param frame An `fmri_frame` or synchronized `fmri_view` with compatible
#'   design annotations and blocks.
#' @return A `compiled_design` for the new observations.
#' @export
apply_design <- function(x, frame) {
  if (!inherits(x, "compiled_design") ||
      !inherits(x$blueprint, "design_blueprint")) {
    stop("`x` must be a compiled_design with a frozen blueprint.", call. = FALSE)
  }
  if (!inherits(frame, c("fmri_frame", "fmri_view"))) {
    stop("`frame` must be an fmri_frame or fmri_view.", call. = FALSE)
  }
  prepared <- .prepare_design_data(
    frame,
    x$spec,
    frozen_mv_terms = x$blueprint$mv_terms
  )
  .assemble_compiled_design(
    frame,
    x$spec,
    prepared,
    frozen_blueprint = x$blueprint
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

#' @rdname model_matrix
#' @export
component_data <- function(x) {
  if (!inherits(x, "compiled_design")) stop("`x` must be a compiled_design.", call. = FALSE)
  x$components
}

#' @rdname model_matrix
#' @export
grouping_term_data <- function(x) {
  if (!inherits(x, "compiled_design")) stop("`x` must be a compiled_design.", call. = FALSE)
  x$grouping_terms
}

#' @rdname model_matrix
#' @export
random_effect_data <- function(x) {
  if (!inherits(x, "compiled_design")) stop("`x` must be a compiled_design.", call. = FALSE)
  x$random_effects
}

#' @rdname model_matrix
#' @export
design_blueprint <- function(x) {
  if (!inherits(x, "compiled_design")) stop("`x` must be a compiled_design.", call. = FALSE)
  x$blueprint
}

#' @rdname model_matrix
#' @export
design_rows <- function(x) {
  if (!inherits(x, "compiled_design")) stop("`x` must be a compiled_design.", call. = FALSE)
  x$row_audit
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
