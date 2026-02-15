#' Combine multiple multidesign objects
#'
#' This function row-binds the observation matrices and design data frames
#' of several multidesign objects. All input multidesigns must share the
#' same column design. Optionally, an identifier column can be added to
#' track the source of each observation.
#'
#' @param ... multidesign objects to combine
#' @param .id optional name of an identifier column added to the design
#'   to indicate the origin multidesign
#'
#' @return A new multidesign object containing the concatenated data
#'
#' @examples
#' X1 <- matrix(rnorm(10*5), 10, 5)
#' X2 <- matrix(rnorm(10*5), 10, 5)
#' md1 <- multidesign(X1, data.frame(cond = rep("A", 10)))
#' md2 <- multidesign(X2, data.frame(cond = rep("B", 10)))
#' combined <- bind_multidesign(md1, md2)
#' combined_id <- bind_multidesign(md1, md2, .id = "source")
#'
#' @export
bind_multidesign <- function(..., .id = NULL) {
  mds_list <- list(...)
  if (length(mds_list) == 1 && is.list(mds_list[[1]])) {
    mds_list <- mds_list[[1]]
  }
  if (length(mds_list) < 1) {
    stop("no multidesign objects supplied")
  }

  # Expand any hyperdesign inputs into their constituent multidesign blocks
  expanded <- list()
  for (item in mds_list) {
    if (inherits(item, "hyperdesign")) {
      expanded <- c(expanded, as.list(item))
    } else {
      expanded <- c(expanded, list(item))
    }
  }
  mds_list <- expanded

  if (!all(sapply(mds_list, inherits, "multidesign"))) {
    stop("all inputs must be multidesign objects")
  }

  base_cd <- mds_list[[1]]$column_design
  for (i in seq_along(mds_list)[-1]) {
    if (!identical(base_cd, mds_list[[i]]$column_design)) {
      stop("column designs must be identical across multidesigns")
    }
  }

  X <- do.call(rbind, lapply(mds_list, function(md) md$x))
  design_list <- lapply(seq_along(mds_list), function(i) {
    des <- mds_list[[i]]$design
    if (!is.null(.id)) des[[.id]] <- i
    des
  })
  design_df <- dplyr::bind_rows(design_list)
  multidesign(X, design_df, base_cd)
}
