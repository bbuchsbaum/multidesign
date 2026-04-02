#' Create Cross-validation Folds from Explicit Row Indices
#'
#' @description
#' Creates cross-validation folds from explicit assessment row indices rather than
#' deriving folds from design variables. This is useful when row holdouts have
#' already been chosen and you want to reuse the package's existing `foldlist`
#' and `cross_validate()` machinery.
#'
#' @param x The dataset to fold over (`multidesign`, `hyperdesign`, or `multiframe`).
#' @param rows Assessment row indices for each fold.
#'   * For `multidesign` and `multiframe`, supply a list of integer vectors, one per fold.
#'   * For `hyperdesign`, supply a list of per-fold block mappings. Each fold is a list whose
#'     elements are integer vectors of row indices, keyed by block name or block position.
#' @param preserve_row_ids Logical; if `TRUE`, carry original source row ids into
#'   fold `analysis` and `assessment` designs via a reserved `.orig_index` column.
#'   Where `held_out` metadata is present, matching `row_ids` are also included.
#' @param ... Additional arguments passed to methods (currently unused).
#'
#' @return A `foldlist` object containing `analysis`, `assessment`, and optional
#'   `held_out` metadata.
#'
#' @examples
#' X <- matrix(rnorm(40), 10, 4)
#' Y <- data.frame(group = rep(c("A", "B"), each = 5))
#' mds <- multidesign(X, Y)
#'
#' folds <- cv_rows(mds, rows = list(1:2, 6:7))
#' folds[[1]]$assessment
#' folds_with_ids <- cv_rows(mds, rows = list(1:2), preserve_row_ids = TRUE)
#' folds_with_ids[[1]]$assessment$design$.orig_index
#'
#' d1 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(run = 1:6))
#' d2 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(run = 1:6))
#' hd <- hyperdesign(list(d1, d2), block_names = c("subj1", "subj2"))
#'
#' hd_folds <- cv_rows(hd, rows = list(
#'   list(subj1 = 1:2, subj2 = 1:2),
#'   list(subj1 = 3:4, subj2 = 3:4)
#' ))
#' hd_folds[[1]]$assessment
#'
#' @seealso
#'   \code{\link{fold_over}} for creating folds from design variables,
#'   \code{\link{cross_validate}} for executing the fit/score loop
#' @export
cv_rows <- function(x, rows, ...) UseMethod("cv_rows")

new_foldlist <- function(extract_fn, fold_ids, foldframe) {
  ret <- deflist::deflist(function(i) extract_fn(fold_ids[[i]]), len = length(fold_ids))
  names(ret) <- paste0("fold_", seq_along(ret))
  class(ret) <- c("foldlist", class(ret))
  attr(ret, "foldframe") <- foldframe
  ret
}

normalize_row_index_vector <- function(indices, n_rows, label, allow_full = FALSE) {
  if (!is.numeric(indices)) {
    stop(label, " must be a numeric or integer vector of row indices.")
  }
  if (length(indices) == 0) {
    stop(label, " must contain at least one row index.")
  }
  if (anyNA(indices) || any(indices != as.integer(indices))) {
    stop(label, " must contain whole-number row indices.")
  }

  indices <- as.integer(indices)

  if (anyDuplicated(indices)) {
    stop(label, " must not contain duplicate row indices.")
  }
  if (any(indices < 1L | indices > n_rows)) {
    stop(label, " contains indices outside the valid range 1:", n_rows, ".")
  }
  if (!allow_full && length(indices) == n_rows) {
    stop(label, " cannot hold out all rows from the dataset.")
  }

  indices
}

normalize_simple_row_folds <- function(rows, n_rows, object_label) {
  if (!is.list(rows) || length(rows) == 0) {
    stop("`rows` must be a non-empty list of assessment row indices for ", object_label, ".")
  }

  lapply(seq_along(rows), function(i) {
    normalize_row_index_vector(
      rows[[i]],
      n_rows = n_rows,
      label = paste0("rows[[", i, "]]")
    )
  })
}

source_row_ids <- function(design_subset) {
  if (".orig_index" %in% names(design_subset)) {
    as.integer(design_subset$.orig_index)
  } else {
    as.integer(design_subset$.index)
  }
}

prepare_design_row_ids <- function(design_subset, preserve_row_ids = FALSE) {
  row_ids <- source_row_ids(design_subset)

  if (preserve_row_ids) {
    design_subset$.orig_index <- row_ids
  } else if (".orig_index" %in% names(design_subset)) {
    design_subset$.orig_index <- NULL
  }

  list(design = design_subset, row_ids = row_ids)
}

slice_multidesign_rows <- function(x, indices, preserve_row_ids = FALSE) {
  prepared <- prepare_design_row_ids(x$design[indices, , drop = FALSE], preserve_row_ids = preserve_row_ids)
  design_subset <- prepared$design
  design_subset$.index <- NULL
  multidesign(x$x[indices, , drop = FALSE], design_subset, x$column_design)
}

exclude_multidesign_rows <- function(x, indices, drop_empty = FALSE, preserve_row_ids = FALSE) {
  keep <- setdiff(seq_len(nrow(x$x)), indices)
  if (length(keep) == 0) {
    if (drop_empty) {
      return(NULL)
    }
    stop("Each fold must leave at least one training row in every affected block.")
  }
  slice_multidesign_rows(x, keep, preserve_row_ids = preserve_row_ids)
}

slice_multiframe_rows <- function(x, indices, preserve_row_ids = FALSE) {
  prepared <- prepare_design_row_ids(x$design[indices, , drop = FALSE], preserve_row_ids = preserve_row_ids)
  design_subset <- prepared$design
  design_subset$.index <- seq_len(nrow(design_subset))
  structure(list(design = design_subset), class = "multiframe")
}

exclude_multiframe_rows <- function(x, indices, preserve_row_ids = FALSE) {
  keep <- setdiff(seq_len(nrow(x$design)), indices)
  if (length(keep) == 0) {
    stop("Each fold must leave at least one training row in the multiframe.")
  }
  slice_multiframe_rows(x, keep, preserve_row_ids = preserve_row_ids)
}

merge_held_out_row_ids <- function(held_out, row_ids, preserve_row_ids = FALSE) {
  if (!preserve_row_ids) {
    return(held_out)
  }

  if (is.null(held_out)) {
    return(list(row_ids = row_ids))
  }
  if (is.list(held_out) && !("row_ids" %in% names(held_out))) {
    held_out$row_ids <- row_ids
  }
  held_out
}

build_multidesign_foldlist <- function(x, foldframe, held_out = NULL, preserve_row_ids = FALSE) {
  fold_ids <- unique(foldframe$.fold)
  if (is.null(held_out)) {
    held_out <- vector("list", length(fold_ids))
  }

  extract <- function(fold_id) {
    fold_row <- foldframe[foldframe$.fold == fold_id, , drop = FALSE]
    if (nrow(fold_row) != 1) {
      stop("Each multidesign fold must define exactly one assessment row set.")
    }

    indices <- fold_row$indices[[1]]
    fold <- list(
      analysis = exclude_multidesign_rows(x, indices, preserve_row_ids = preserve_row_ids),
      assessment = slice_multidesign_rows(x, indices, preserve_row_ids = preserve_row_ids)
    )

    fold_index <- match(fold_id, fold_ids)
    held_out_info <- held_out[[fold_index]]
    row_ids <- source_row_ids(x$design[indices, , drop = FALSE])
    held_out_info <- merge_held_out_row_ids(held_out_info, row_ids, preserve_row_ids = preserve_row_ids)
    if (!is.null(held_out_info)) {
      fold$held_out <- held_out_info
    }

    fold
  }

  new_foldlist(extract, fold_ids, foldframe)
}

build_multiframe_foldlist <- function(x, foldframe, held_out = NULL, preserve_row_ids = FALSE) {
  fold_ids <- unique(foldframe$.fold)
  if (is.null(held_out)) {
    held_out <- vector("list", length(fold_ids))
  }

  extract <- function(fold_id) {
    fold_row <- foldframe[foldframe$.fold == fold_id, , drop = FALSE]
    if (nrow(fold_row) != 1) {
      stop("Each multiframe fold must define exactly one assessment row set.")
    }

    indices <- fold_row$indices[[1]]
    fold <- list(
      analysis = exclude_multiframe_rows(x, indices, preserve_row_ids = preserve_row_ids),
      assessment = slice_multiframe_rows(x, indices, preserve_row_ids = preserve_row_ids)
    )

    fold_index <- match(fold_id, fold_ids)
    held_out_info <- held_out[[fold_index]]
    row_ids <- source_row_ids(x$design[indices, , drop = FALSE])
    held_out_info <- merge_held_out_row_ids(held_out_info, row_ids, preserve_row_ids = preserve_row_ids)
    if (!is.null(held_out_info)) {
      fold$held_out <- held_out_info
    }

    fold
  }

  new_foldlist(extract, fold_ids, foldframe)
}

resolve_hyperdesign_block_positions <- function(x, spec, fold_label) {
  block_names <- names(x)
  if (is.null(block_names)) {
    block_names <- paste0("block_", seq_along(x))
  }

  spec_names <- names(spec)
  if (is.null(spec_names) || any(spec_names == "")) {
    positions <- seq_along(spec)
  } else {
    positions <- vapply(seq_along(spec), function(i) {
      nm <- spec_names[[i]]
      if (!is.na(suppressWarnings(as.integer(nm))) && nm %in% as.character(seq_along(x))) {
        as.integer(nm)
      } else {
        match(nm, block_names)
      }
    }, integer(1))

    if (anyNA(positions)) {
      stop(fold_label, " refers to unknown hyperdesign blocks.")
    }
  }

  if (any(positions < 1L | positions > length(x))) {
    stop(fold_label, " refers to block indices outside the valid range 1:", length(x), ".")
  }
  if (anyDuplicated(positions)) {
    stop(fold_label, " must not reference the same block more than once.")
  }

  list(positions = positions, block_names = block_names)
}

normalize_hyperdesign_row_folds <- function(x, rows, allow_full_block_holdout = FALSE) {
  if (!is.list(rows) || length(rows) == 0) {
    stop("`rows` must be a non-empty list of per-fold block row mappings for hyperdesign objects.")
  }

  block_names <- names(x)
  if (is.null(block_names)) {
    block_names <- paste0("block_", seq_along(x))
  }

  foldframe_parts <- vector("list", length(rows))
  held_out <- vector("list", length(rows))

  for (i in seq_along(rows)) {
    fold_label <- paste0("rows[[", i, "]]")
    spec <- rows[[i]]

    if (!is.list(spec) || length(spec) == 0) {
      stop(fold_label, " must be a non-empty list of block-specific row indices.")
    }

    resolved <- resolve_hyperdesign_block_positions(x, spec, fold_label)
    positions <- resolved$positions

    normalized_indices <- lapply(seq_along(spec), function(j) {
      block_pos <- positions[[j]]
      normalize_row_index_vector(
        spec[[j]],
        n_rows = nrow(x[[block_pos]]$x),
        label = paste0(fold_label, "[[", j, "]]"),
        allow_full = allow_full_block_holdout
      )
    })

    held_out[[i]] <- stats::setNames(normalized_indices, block_names[positions])
    foldframe_parts[[i]] <- tibble::tibble(
      .fold = i,
      .block = positions,
      indices = normalized_indices,
      .splitvar = paste(block_names[positions], collapse = "+")
    )
  }

  list(
    foldframe = dplyr::bind_rows(foldframe_parts),
    held_out = held_out
  )
}

build_hyperdesign_foldlist <- function(x,
                                       foldframe,
                                       held_out = NULL,
                                       assessment_mode = c("multidesign", "hyperdesign"),
                                       drop_empty_analysis_blocks = FALSE,
                                       preserve_row_ids = FALSE) {
  assessment_mode <- match.arg(assessment_mode)
  fold_ids <- unique(foldframe$.fold)
  block_names <- names(x)
  if (is.null(block_names)) {
    block_names <- paste0("block_", seq_along(x))
  }

  if (is.null(held_out)) {
    held_out <- vector("list", length(fold_ids))
  }

  extract <- function(fold_id) {
    fold_rows <- foldframe[foldframe$.fold == fold_id, , drop = FALSE]
    assessment_blocks <- lapply(seq_len(nrow(fold_rows)), function(i) {
      block_pos <- fold_rows$.block[[i]]
      slice_multidesign_rows(x[[block_pos]], fold_rows$indices[[i]], preserve_row_ids = preserve_row_ids)
    })
    assessment_names <- block_names[fold_rows$.block]

    analysis_blocks <- lapply(seq_along(x), function(i) {
      match_idx <- which(fold_rows$.block == i)
      if (length(match_idx) == 0) {
        x[[i]]
      } else {
        exclude_multidesign_rows(
          x[[i]],
          fold_rows$indices[[match_idx[[1]]]],
          drop_empty = drop_empty_analysis_blocks,
          preserve_row_ids = preserve_row_ids
        )
      }
    })

    if (drop_empty_analysis_blocks) {
      keep <- !vapply(analysis_blocks, is.null, logical(1))
      analysis_blocks <- analysis_blocks[keep]
      analysis_names <- block_names[keep]
    } else {
      analysis_names <- block_names
    }

    if (length(analysis_blocks) == 0) {
      stop("Each fold must leave at least one analysis block in the hyperdesign.")
    }

    fold <- list(
      analysis = hyperdesign(analysis_blocks, block_names = analysis_names),
      assessment = if (assessment_mode == "multidesign") {
        if (length(assessment_blocks) != 1) {
          stop("Single-block assessment mode requires exactly one assessment block per fold.")
        }
        assessment_blocks[[1]]
      } else {
        hyperdesign(assessment_blocks, block_names = assessment_names)
      }
    )

    fold_index <- match(fold_id, fold_ids)
    held_out_info <- held_out[[fold_index]]
    row_ids <- stats::setNames(lapply(seq_len(nrow(fold_rows)), function(i) {
      block_pos <- fold_rows$.block[[i]]
      source_row_ids(x[[block_pos]]$design[fold_rows$indices[[i]], , drop = FALSE])
    }), assessment_names)
    held_out_info <- merge_held_out_row_ids(held_out_info, row_ids, preserve_row_ids = preserve_row_ids)
    if (!is.null(held_out_info)) {
      fold$held_out <- held_out_info
    }

    fold
  }

  new_foldlist(extract, fold_ids, foldframe)
}

#' @rdname cv_rows
#' @export
cv_rows.multidesign <- function(x, rows, preserve_row_ids = FALSE, ...) {
  rows <- normalize_simple_row_folds(rows, nrow(x$x), "multidesign")
  foldframe <- tibble::tibble(
    indices = rows,
    .splitvar = paste0("rows_", seq_along(rows)),
    .fold = seq_along(rows)
  )
  held_out <- lapply(rows, function(indices) list(rows = indices))
  build_multidesign_foldlist(x, foldframe, held_out = held_out, preserve_row_ids = preserve_row_ids)
}

#' @rdname cv_rows
#' @export
cv_rows.multiframe <- function(x, rows, preserve_row_ids = FALSE, ...) {
  rows <- normalize_simple_row_folds(rows, nrow(x$design), "multiframe")
  foldframe <- tibble::tibble(
    indices = rows,
    .splitvar = paste0("rows_", seq_along(rows)),
    .fold = seq_along(rows)
  )
  held_out <- lapply(rows, function(indices) list(rows = indices))
  build_multiframe_foldlist(x, foldframe, held_out = held_out, preserve_row_ids = preserve_row_ids)
}

#' @rdname cv_rows
#' @export
cv_rows.hyperdesign <- function(x, rows, preserve_row_ids = FALSE, ...) {
  normalized <- normalize_hyperdesign_row_folds(x, rows)
  build_hyperdesign_foldlist(
    x,
    normalized$foldframe,
    held_out = normalized$held_out,
    assessment_mode = "hyperdesign",
    preserve_row_ids = preserve_row_ids
  )
}
