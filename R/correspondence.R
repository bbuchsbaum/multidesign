#' Query Row Correspondence in a Hyperdesign
#'
#' A correspondence object describes how local rows in each hyperdesign block
#' map to a global entity universe. Correspondence is opt-in: construct the
#' hyperdesign with either `id` or `positional = TRUE`.
#'
#' Explicit identifiers are compared through a common character representation,
#' so integer `1` and character `"1"` match. `global_ids` follows first-seen
#' block and row order, and each integer `row_map` indexes that vector. Missing
#' values in the data matrix are ordinary data and are never interpreted as a
#' correspondence mask.
#'
#' @param x A hyperdesign object.
#' @param sort_ids Logical; if `TRUE`, order normalized global identifiers with
#'   base R radix ordering. The default preserves first-seen order.
#' @param ... Additional arguments passed to methods.
#'
#' @return For a contracted hyperdesign, a `md_correspondence` list containing
#'   the declared key, assumption, global identifiers, local identifiers,
#'   integer row maps, incidence matrix, pairwise overlap counts, and block
#'   connectivity flag.
#' @examples
#' a <- multidesign(matrix(1:6, ncol = 2), data.frame(id = c("B", "A", "C")))
#' b <- multidesign(matrix(7:10, ncol = 2), data.frame(id = c("A", "D")))
#' hd <- hyperdesign(list(a = a, b = b), id = "id", space = "common")
#' correspondence(hd)$global_ids
#' correspondence(hd, sort_ids = TRUE)$global_ids
#' @export
correspondence <- function(x, ...) UseMethod("correspondence")

#' @rdname correspondence
#' @method correspondence hyperdesign
#' @export
correspondence.hyperdesign <- function(x, sort_ids = FALSE, ...) {
  if (!is.logical(sort_ids) || length(sort_ids) != 1L || is.na(sort_ids)) {
    stop("`sort_ids` must be TRUE or FALSE.", call. = FALSE)
  }
  id <- entity_id(x)
  space <- column_space(x)
  assumption <- attr(x, "correspondence_assumption", exact = TRUE)

  if (!is.null(id)) {
    assumption <- "explicit_ids"
  }
  if (is.null(assumption)) {
    stop(
      "`x` has no row-correspondence contract; reconstruct it with `id =` or `positional = TRUE`.",
      call. = FALSE
    )
  }
  if (!assumption %in% c("explicit_ids", "positional")) {
    stop("`x` has an invalid correspondence assumption.", call. = FALSE)
  }

  contract <- .validate_hyperdesign_contract(
    x,
    id = id,
    space = space,
    positional = identical(assumption, "positional")
  )

  if (identical(assumption, "explicit_ids")) {
    if (is.null(contract$id)) {
      stop("Explicit correspondence requires a valid entity ID column.", call. = FALSE)
    }
    ids <- contract$ids
  } else {
    global <- as.character(seq_len(nrow(x[[1L]]$x)))
    ids <- lapply(x, function(block) global)
    names(ids) <- names(x)
  }

  global_ids <- unique(unlist(ids, use.names = FALSE))
  if (isTRUE(sort_ids)) {
    global_ids <- sort(global_ids, method = "radix")
  }
  row_map <- lapply(ids, match, table = global_ids)
  names(row_map) <- names(x)

  incidence <- matrix(
    FALSE,
    nrow = length(global_ids),
    ncol = length(x),
    dimnames = list(global_ids, names(x))
  )
  for (i in seq_along(row_map)) {
    incidence[row_map[[i]], i] <- TRUE
  }

  pair_n <- crossprod(incidence)
  storage.mode(pair_n) <- "integer"
  connected <- .correspondence_connected(pair_n > 0L)

  out <- list(
    id = contract$id,
    assumption = assumption,
    global_ids = global_ids,
    ids = ids,
    row_map = row_map,
    n_global = length(global_ids),
    n_observed = stats::setNames(
      vapply(ids, length, integer(1)),
      names(x)
    ),
    overlap = list(
      incidence = incidence,
      pair_n = pair_n,
      connected = connected
    )
  )
  class(out) <- c("md_correspondence", "list")
  out
}

#' Test Whether a Hyperdesign Declares Row Correspondence
#'
#' @param x A hyperdesign object.
#' @return A length-one logical value.
#' @export
has_correspondence <- function(x) {
  .check_hyperdesign(x)
  !is.null(entity_id(x)) ||
    identical(attr(x, "correspondence_assumption", exact = TRUE), "positional")
}

#' Get the Entity ID Column of a Hyperdesign
#'
#' @param x A hyperdesign object.
#' @return The design-column name, or `NULL` when unset or positional.
#' @export
entity_id <- function(x) {
  .check_hyperdesign(x)
  attr(x, "entity_id", exact = TRUE)
}

#' Get the Column-Space Contract of a Hyperdesign
#'
#' @param x A hyperdesign object.
#' @return `"common"`, `"block"`, or `NULL` when unspecified.
#' @export
column_space <- function(x) {
  .check_hyperdesign(x)
  attr(x, "column_space", exact = TRUE)
}

.check_hyperdesign <- function(x) {
  if (!inherits(x, "hyperdesign")) {
    stop("`x` must inherit from `hyperdesign`.", call. = FALSE)
  }
  invisible(x)
}

.validate_hyperdesign_contract <- function(x, id = NULL, space = NULL,
                                           positional = FALSE) {
  if (!is.logical(positional) || length(positional) != 1L || is.na(positional)) {
    stop("`positional` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(id)) {
    if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
      stop("`id` must be NULL or one nonempty design-column name.", call. = FALSE)
    }
    if (id %in% c(".index", ".orig_index")) {
      stop("`id` cannot use reserved row-index columns `.index` or `.orig_index`.", call. = FALSE)
    }
    if (isTRUE(positional)) {
      stop("`id` and `positional = TRUE` are mutually exclusive.", call. = FALSE)
    }
  }

  if (!is.null(space)) {
    if (!is.character(space) || length(space) != 1L || is.na(space) ||
        !space %in% c("common", "block")) {
      stop("`space` must be NULL, \"common\", or \"block\".", call. = FALSE)
    }
  }

  contracted <- !is.null(id) || !is.null(space) || isTRUE(positional)
  if (contracted) {
    block_names <- names(x)
    if (is.null(block_names) || anyNA(block_names) || any(!nzchar(block_names)) ||
        anyDuplicated(block_names)) {
      stop("Contracted hyperdesign blocks must have unique, nonempty names.", call. = FALSE)
    }
  }

  ids <- NULL
  if (!is.null(id)) {
    ids <- .normalize_hyperdesign_ids(x, id)
  }

  if (isTRUE(positional)) {
    nrows <- vapply(x, function(block) nrow(block$x), integer(1))
    if (length(unique(nrows)) != 1L) {
      stop("`positional = TRUE` requires equal row counts in every block.", call. = FALSE)
    }
  }

  if (identical(space, "common")) {
    ncols <- vapply(x, function(block) ncol(block$x), integer(1))
    if (length(unique(ncols)) != 1L) {
      stop("`space = \"common\"` requires equal column counts in every block.", call. = FALSE)
    }
    base_column_design <- x[[1L]]$column_design
    if (length(x) > 1L) {
      for (i in 2:length(x)) {
        if (!identical(base_column_design, x[[i]]$column_design)) {
          stop(
            "`space = \"common\"` requires identical column designs across blocks; block ",
            i,
            " differs.",
            call. = FALSE
          )
        }
      }
    }
  }

  list(id = id, space = space, positional = positional, ids = ids)
}

.normalize_hyperdesign_ids <- function(x, id) {
  out <- lapply(seq_along(x), function(i) {
    design <- x[[i]]$design
    if (!id %in% names(design)) {
      stop(
        "Entity ID column `", id, "` is missing from block `", names(x)[[i]], "`.",
        call. = FALSE
      )
    }
    values <- design[[id]]
    if (!is.atomic(values) || !is.null(dim(values))) {
      stop(
        "Entity ID column `", id, "` in block `", names(x)[[i]], "` must be an atomic vector.",
        call. = FALSE
      )
    }
    if (length(values) != nrow(x[[i]]$x)) {
      stop(
        "Entity ID column `", id, "` in block `", names(x)[[i]], "` must have one value per row.",
        call. = FALSE
      )
    }
    if (anyNA(values)) {
      stop(
        "Entity ID column `", id, "` in block `", names(x)[[i]], "` cannot contain NA.",
        call. = FALSE
      )
    }
    normalized <- as.character(values)
    if (anyNA(normalized)) {
      stop(
        "Entity ID column `", id, "` in block `", names(x)[[i]], "` cannot normalize to NA.",
        call. = FALSE
      )
    }
    if (anyDuplicated(normalized)) {
      stop(
        "Entity ID column `", id, "` contains duplicate keys in block `", names(x)[[i]], "`.",
        call. = FALSE
      )
    }
    unname(normalized)
  })
  names(out) <- names(x)
  out
}

.correspondence_connected <- function(adjacency) {
  k <- nrow(adjacency)
  if (k <= 1L) {
    return(TRUE)
  }
  diag(adjacency) <- TRUE
  seen <- rep(FALSE, k)
  queue <- 1L
  while (length(queue)) {
    current <- queue[[1L]]
    queue <- queue[-1L]
    if (seen[[current]]) {
      next
    }
    seen[[current]] <- TRUE
    neighbours <- which(adjacency[current, ] & !seen)
    queue <- unique(c(queue, neighbours))
  }
  all(seen)
}
