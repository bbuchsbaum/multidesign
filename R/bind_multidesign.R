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
#' @export
bind_multidesign <- function(..., .id = NULL) {
  mds_list <- list(...)
  if (length(mds_list) == 1 && is.list(mds_list[[1]])) {
    mds_list <- mds_list[[1]]
  }
  if (length(mds_list) < 1) {
    stop("no multidesign objects supplied")
  }
  if (!all(sapply(mds_list, inherits, "multidesign"))) {
    stop("all inputs must be multidesign objects")
  }

  base_cd <- mds_list[[1]]$column_design
  for (md in mds_list[-1]) {
    if (!identical(base_cd, md$column_design)) {
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
