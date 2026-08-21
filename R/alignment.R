#' Align Hyperdesign Blocks by Their Entity Contract
#'
#' Places every common-space block on the global entity universe described by
#' [correspondence()]. The result keeps row presence and cell observation as
#' separate logical arrays; the fill value is presentation only.
#'
#' @param x A hyperdesign with declared correspondence and `space = "common"`.
#' @param fill Length-one numeric value used to display absent rows and
#'   explicitly unobserved cells.
#' @param sort_ids Logical; if `TRUE`, use sorted global identifiers rather than
#'   first-seen order.
#' @param ... Additional arguments passed to methods.
#' @return An `aligned_hyperdesign` list containing `x`, `cells`, `observed`,
#'   `correspondence`, `column_design`, and `block_names`.
#' @examples
#' a <- multidesign(
#'   matrix(c(NA, 0, 1, 0), ncol = 2, byrow = TRUE),
#'   data.frame(id = c("A", "B")),
#'   cells = matrix(c(TRUE, TRUE, TRUE, FALSE), ncol = 2, byrow = TRUE)
#' )
#' b <- multidesign(
#'   matrix(c(0, 1, 1, 1), ncol = 2, byrow = TRUE),
#'   data.frame(id = c("B", "C"))
#' )
#' hd <- hyperdesign(list(a = a, b = b), id = "id", space = "common")
#' aligned <- align_by_id(hd)
#' aligned$observed
#' aligned$cells[, , "a"]
#' aligned$x[, , "a"]
#' @export
align_by_id <- function(x, ...) UseMethod("align_by_id")

#' @rdname align_by_id
#' @method align_by_id hyperdesign
#' @export
align_by_id.hyperdesign <- function(x, fill = NA_real_, sort_ids = FALSE, ...) {
  if (!has_correspondence(x)) {
    stop(
      "`align_by_id()` requires declared correspondence via `id =` or `positional = TRUE`.",
      call. = FALSE
    )
  }
  if (!identical(column_space(x), "common")) {
    stop(
      "`align_by_id()` requires `space = \"common\"`; block or unspecified column spaces cannot be aligned.",
      call. = FALSE
    )
  }
  if (!is.numeric(fill) || length(fill) != 1L) {
    stop("`fill` must be one numeric value.", call. = FALSE)
  }

  cr <- correspondence(x, sort_ids = sort_ids)
  n_entities <- cr$n_global
  n_variables <- ncol(x[[1L]]$x)
  n_blocks <- length(x)
  block_names <- names(x)

  variable_names <- colnames(x[[1L]]$x)
  if (is.null(variable_names)) {
    variable_names <- as.character(seq_len(n_variables))
  }
  array_dimnames <- list(
    entity = cr$global_ids,
    variable = variable_names,
    block = block_names
  )
  aligned_x <- array(
    as.numeric(fill),
    dim = c(n_entities, n_variables, n_blocks),
    dimnames = array_dimnames
  )
  aligned_cells <- array(
    FALSE,
    dim = c(n_entities, n_variables, n_blocks),
    dimnames = array_dimnames
  )

  for (block_index in seq_along(x)) {
    block <- x[[block_index]]
    row_map <- cr$row_map[[block_index]]
    mask <- cell_mask(block)
    if (is.null(mask)) {
      mask <- matrix(TRUE, nrow = nrow(block$x), ncol = ncol(block$x))
    }
    presented <- matrix(
      as.numeric(fill),
      nrow = nrow(block$x),
      ncol = ncol(block$x)
    )
    presented[mask] <- block$x[mask]

    aligned_x[row_map, , block_index] <- presented
    aligned_cells[row_map, , block_index] <- mask
  }

  structure(
    list(
      x = aligned_x,
      cells = aligned_cells,
      observed = cr$overlap$incidence,
      correspondence = cr,
      column_design = x[[1L]]$column_design,
      block_names = block_names
    ),
    class = c("aligned_hyperdesign", "list")
  )
}

#' Print an Aligned Hyperdesign
#'
#' @param x An aligned_hyperdesign object.
#' @param ... Additional arguments (unused).
#' @return `x`, invisibly.
#' @method print aligned_hyperdesign
#' @export
print.aligned_hyperdesign <- function(x, ...) {
  dimensions <- dim(x$x)
  cat("\nAligned hyperdesign\n")
  cat("  ", dimensions[[1L]], " entities x ", dimensions[[2L]],
      " variables x ", dimensions[[3L]], " blocks\n", sep = "")
  cat("  Row observations: ", sum(x$observed), " of ", length(x$observed), "\n", sep = "")
  cat("  Cell observations: ", sum(x$cells), " of ", length(x$cells), "\n", sep = "")
  invisible(x)
}
