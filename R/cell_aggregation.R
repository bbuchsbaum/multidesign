.resolve_cell_aggregate <- function(aggregate) {
  if (is.character(aggregate) && length(aggregate) == 1L &&
      !is.na(aggregate) && identical(aggregate, "mean")) {
    return(base::mean)
  }
  if (is.function(aggregate)) {
    return(aggregate)
  }
  stop(
    "`aggregate` must be NULL, \"mean\", or a scalar-returning function.",
    call. = FALSE
  )
}

.aggregate_matrix_groups <- function(data, cells, row_groups, aggregate,
                                     labels, context,
                                     preserve_singletons = FALSE) {
  aggregate_fun <- .resolve_cell_aggregate(aggregate)
  n_groups <- length(row_groups)
  n_variables <- ncol(data)
  explicit_mask <- !is.null(cells)
  observed <- if (explicit_mask) {
    cells
  } else {
    matrix(TRUE, nrow = nrow(data), ncol = ncol(data))
  }

  out_x <- matrix(
    NA_real_,
    nrow = n_groups,
    ncol = n_variables,
    dimnames = list(labels, colnames(data))
  )
  out_cells <- matrix(
    FALSE,
    nrow = n_groups,
    ncol = n_variables,
    dimnames = list(labels, colnames(data))
  )

  for (group_index in seq_along(row_groups)) {
    rows <- as.integer(row_groups[[group_index]])
    if (length(rows) == 0L) {
      stop("Aggregation groups must contain at least one row.", call. = FALSE)
    }

    if (isTRUE(preserve_singletons) && length(rows) == 1L) {
      out_x[group_index, ] <- data[rows, , drop = TRUE]
      out_cells[group_index, ] <- observed[rows, , drop = TRUE]
      next
    }

    for (coordinate in seq_len(n_variables)) {
      coordinate_observed <- observed[rows, coordinate]
      if (!any(coordinate_observed)) {
        next
      }

      values <- data[rows[coordinate_observed], coordinate]
      value <- tryCatch(
        aggregate_fun(values),
        error = function(error) {
          stop(
            "Aggregation failed for ", context, " `", labels[[group_index]],
            "`, coordinate ", coordinate, ": ", conditionMessage(error),
            call. = FALSE
          )
        }
      )
      if (!is.numeric(value) || length(value) != 1L) {
        stop(
          "Aggregation for ", context, " `", labels[[group_index]],
          "`, coordinate ", coordinate,
          " must return one numeric value.",
          call. = FALSE
        )
      }
      out_x[group_index, coordinate] <- as.numeric(value)
      out_cells[group_index, coordinate] <- TRUE
    }
  }

  list(x = out_x, cells = if (explicit_mask) out_cells else NULL)
}

.design_values_identical <- function(values) {
  if (length(values) <= 1L) {
    return(TRUE)
  }
  reference <- values[[1L]]
  all(vapply(values[-1L], identical, logical(1), y = reference))
}

.aggregate_duplicate_entity_block <- function(block, id, aggregate,
                                               block_name) {
  normalized_ids <- as.character(block$design[[id]])
  group_labels <- unique(normalized_ids)
  row_groups <- lapply(group_labels, function(label) which(normalized_ids == label))

  if (!any(vapply(row_groups, length, integer(1)) > 1L)) {
    return(block)
  }

  design_columns <- setdiff(names(block$design), c(".index", ".orig_index"))
  metadata_columns <- setdiff(design_columns, id)
  for (group_index in seq_along(row_groups)) {
    rows <- row_groups[[group_index]]
    if (length(rows) <= 1L) {
      next
    }
    conflicting <- metadata_columns[!vapply(
      metadata_columns,
      function(column) .design_values_identical(block$design[[column]][rows]),
      logical(1)
    )]
    if (length(conflicting) > 0L) {
      stop(
        "Cannot aggregate entity `", group_labels[[group_index]],
        "` in block `", block_name,
        "`: non-key design fields disagree: ",
        paste(conflicting, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  aggregated <- .aggregate_matrix_groups(
    block$x,
    cells = cell_mask(block),
    row_groups = row_groups,
    aggregate = aggregate,
    labels = group_labels,
    context = paste0("block ", block_name, ", entity"),
    preserve_singletons = TRUE
  )
  first_rows <- vapply(row_groups, `[[`, integer(1), 1L)
  design_out <- block$design[first_rows, design_columns, drop = FALSE]

  multidesign(
    aggregated$x,
    design_out,
    block$column_design,
    cells = aggregated$cells
  )
}
