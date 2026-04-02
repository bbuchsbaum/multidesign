#' @importFrom magrittr %>%
#' @importFrom utils head tail
#' @importFrom dplyr filter ungroup pull all_of
#' @title Create a Multidesign Object from Matrix Data and Design Information
#'
#' @description
#' Creates a multidesign object that combines experimental data (as a matrix) with design information
#' (as a data frame) and optional column metadata. This structure is particularly useful for
#' experimental designs where observations have multiple associated factors and variables may have metadata.
#'
#' @details
#' A multidesign object consists of three main components:
#' * A data matrix where rows represent observations and columns represent variables
#' * A design data frame containing experimental factors and conditions for each observation
#' * Optional column metadata describing properties of each variable
#'
#' The object maintains the relationship between these components while providing methods
#' for manipulation, subsetting, and analysis.
#'
#' @param x A numeric matrix where rows are observations and columns are variables
#' @param y A data frame containing design variables for each observation (must have same number of rows as x)
#' @param column_design Optional data frame containing metadata for columns in x (must have same number of rows as ncol(x))
#'
#' @return A multidesign object with components:
#'   \item{x}{The input data matrix}
#'   \item{design}{A tibble containing design variables with an added .index column}
#'   \item{column_design}{A tibble containing column metadata (empty if not provided)}
#'
#' @examples
#' # Create example data matrix
#' X <- matrix(rnorm(20*100), 20, 100)
#'
#' # Create design information
#' Y <- tibble::tibble(
#'   condition = rep(c("control", "treatment"), each=10),
#'   subject = rep(1:5, times=4)
#' )
#'
#' # Create column metadata
#' col_info <- data.frame(
#'   roi = paste0("region_", 1:100),
#'   hemisphere = rep(c("left", "right"), 50)
#' )
#'
#' # Create multidesign object
#' mds <- multidesign(X, Y, col_info)
#'
#' @family multidesign functions
#' @seealso
#'   \code{\link{reduce.multidesign}} for dimensionality reduction,
#'   \code{\link{split.multidesign}} for splitting by design variables
#' @importFrom dplyr mutate rowwise n
#' @export
#' @rdname multidesign
multidesign.matrix <- function(x, y, column_design=NULL, ...) {
  chk::chk_equal(nrow(x), nrow(y))
  chk::chk_s3_class(y, "data.frame")
  y <- tibble::as_tibble(y)
  y$.index <- seq_len(nrow(x))

  if (!is.null(column_design)) {
    chk::chk_equal(ncol(x), nrow(column_design))
    chk::chk_s3_class(column_design, "data.frame")
    column_design <- tibble::as_tibble(column_design)
  } else {
    column_design <- tibble::tibble(.index = seq_len(ncol(x)))
  }

  structure(list(
    x=x,
    design=y,
    column_design=column_design
  ),
  class="multidesign")
}

#' Reduce Dimensionality of a Multidesign Object
#'
#' @description
#' Performs dimensionality reduction on the data matrix of a multidesign object while preserving
#' the design structure. By default uses PCA, but supports any reduction method that returns
#' a projector object.
#'
#' @param x A multidesign object
#' @param nc Number of components to retain in the reduction
#' @param ... Additional arguments passed to rfun
#' @param rfun Function to perform dimensionality reduction, must return a projector object
#'
#' @return A reduced_multidesign object containing:
#'   \item{x}{The reduced data matrix}
#'   \item{design}{Original design information}
#'   \item{projector}{The projection object used for reduction}
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100*20), 100, 20)
#' Y <- tibble::tibble(condition = rep(letters[1:4], each=25))
#' mds <- multidesign(X, Y)
#'
#' # Reduce to 5 components using PCA
#' reduced_mds <- reduce(mds, nc=5)
#' }
#'
#' @family multidesign functions
#' @seealso \code{\link{multidesign}}
#' @export
reduce.multidesign <- function(x, nc=2, ..., rfun=function(x) multivarious::pca(x$x, ncomp=nc,...)) {
  projector <- rfun(x)
  chk::chk_s3_class(projector, "projector")
  rx <- multivarious::project(projector, x$x)
  structure(list(
    x=rx,
    design=x$design,
    column_design=x$column_design,
    projector=projector
  ),
  class=c("reduced_multidesign", "multidesign"))
}

#' Subset a Multidesign Object
#'
#' @description
#' Creates a new multidesign object containing only observations that meet specified
#' conditions based on design variables.
#'
#' @param x A multidesign object
#' @param fexpr An expression used to filter observations based on design variables
#' @param ... Additional arguments (not used)
#'
#' @return A new multidesign object containing only matching observations, or NULL if no matches
#'
#' @examples
#' X <- matrix(rnorm(100*20), 100, 20)
#' Y <- tibble::tibble(
#'   condition = rep(c("A", "B"), each=50),
#'   block = rep(1:2, times=50)
#' )
#' mds <- multidesign(X, Y)
#'
#' # Subset to condition A only
#' mds_A <- subset(mds, condition == "A")
#'
#' # Subset to block 1 and condition A
#' mds_A1 <- subset(mds, condition == "A" & block == 1)
#'
#' @family multidesign functions
#' @seealso \code{\link{split.multidesign}}
#' @export
subset.multidesign <- function(x, fexpr, ...) {
  des2 <- filter(x$design, !!rlang::enquo(fexpr))
  if (nrow(des2) == 0) {
    NULL
  } else {
    ind <- des2$.index
    des2$.index <- NULL
    multidesign(x$x[ind, , drop=FALSE], des2, x$column_design)
  }
}

#' Split a Multidesign Object by Design Variables
#'
#' @description
#' Splits a multidesign object into multiple objects based on combinations of one or more
#' design variables.
#'
#' @param x A multidesign object
#' @param f Unquoted name of the first variable to split by
#' @param drop Logical; unused, present for compatibility with the generic
#' @param ... Additional unquoted names of variables to split by
#'
#' @return A list of multidesign objects, one for each combination of splitting variables
#'
#' @examples
#' X <- matrix(rnorm(100*20), 100, 20)
#' Y <- tibble::tibble(
#'   condition = rep(c("A", "B"), each=50),
#'   block = rep(1:2, times=50)
#' )
#' mds <- multidesign(X, Y)
#'
#' # Split by condition
#' split_by_cond <- split(mds, condition)
#'
#' # Split by condition and block
#' split_by_both <- split(mds, condition, block)
#'
#' @family multidesign functions
#' @seealso \code{\link{subset.multidesign}}
#' @export
split.multidesign <- function(x, f, drop=FALSE, ...) {
  nest.by <- rlang::enquos(f, ...)
  ret <- x$design %>% nest_by(!!!nest.by, .keep=TRUE)
  xl <- ret$data %>% purrr::map(~{
    ind <- .x$.index
    x$x[ind, , drop=FALSE]
  })
  lapply(seq_len(nrow(ret)), function(i) {
    d <- ret$data[[i]]
    d$.index <- NULL
    multidesign.matrix(xl[[i]], d, x$column_design)
  })
}

#' Get Split Indices for a Multidesign Object
#'
#' @description
#' Returns indices for splitting a multidesign object based on combinations of design variables,
#' without actually splitting the data.
#'
#' @param x A multidesign object
#' @param ... Unquoted names of variables to split by
#' @param collapse Logical; whether to collapse the resulting indices
#'
#' @return A tibble containing:
#'   \item{group variables}{The splitting variables and their values}
#'   \item{indices}{List column of indices for each group}
#'   \item{.splitvar}{Combined string of grouping values}
#'
#' @examples
#' X <- matrix(rnorm(40*10), 40, 10)
#' Y <- data.frame(condition = rep(c("A","B"), each=20),
#'                 block = rep(1:4, each=10))
#' mds <- multidesign(X, Y)
#' idx <- split_indices(mds, condition)
#'
#' @family multidesign functions
#' @seealso \code{\link{split.multidesign}}
#' @export
#' @importFrom tidyr unite
#' @importFrom dplyr across
#' @importFrom tidyselect where any_of
split_indices.multidesign <- function(x, ..., collapse=FALSE) {
  nest.by <- rlang::quos(...)

  # Convert numeric variables to factors for proper grouping, excluding row index metadata
  design_copy <- x$design %>%
    mutate(across(where(is.numeric) & !any_of(c(".index", ".orig_index")), as.factor))

  ret <- design_copy %>% nest_by(!!!nest.by, .keep=TRUE)
  xl <- ret$data %>% purrr::map(~ .x$.index)

  # Get group variable names
  group_vars <- colnames(ret %>% select(dplyr::group_vars(ret)))

  # Create .splitvar by uniting all group variables with underscore
  selret <- ret %>% select(dplyr::group_vars(ret))
  out <- selret %>%
    ungroup() %>%
    mutate(
      indices = xl,
      .splitvar = if (length(group_vars) > 0) {
        unite(selret, "temp", all_of(group_vars), sep="_") %>% pull(temp)
      } else {
        rep("all", nrow(selret))
      }
    )

  out
}

#' Summarize a Multidesign Object by Design Variables
#'
#' @description
#' Computes summaries of the data matrix grouped by combinations of design variables.
#'
#' @param x A multidesign object
#' @param ... Unquoted names of variables to group by
#' @param sfun Summary function to apply (default is colMeans)
#' @param extract_data Logical; whether to extract raw data instead of computing summary
#'
#' @return A new multidesign object containing:
#'   \item{x}{Matrix of summary statistics}
#'   \item{design}{Design information for each summary}
#'   \item{column_design}{Original column metadata}
#'
#' @examples
#' X <- matrix(rnorm(100*20), 100, 20)
#' Y <- tibble::tibble(
#'   condition = rep(c("A", "B"), each=50),
#'   block = rep(1:2, times=50)
#' )
#' mds <- multidesign(X, Y)
#'
#' # Get means by condition
#' means_by_cond <- summarize_by(mds, condition)
#'
#' # Get means by condition and block
#' means_by_both <- summarize_by(mds, condition, block)
#'
#' @family multidesign functions
#' @seealso \code{\link{split.multidesign}}
#' @export
summarize_by.multidesign <- function(x, ..., sfun=colMeans, extract_data=FALSE) {
  #nested <- split(x, ...)
  nest.by <- rlang::quos(...)
  ret <- x$design %>% nest_by(!!!nest.by)

  dsum <- do.call(rbind, ret$data %>% purrr::map( ~ sfun(x$x[.x[[".index"]], , drop=FALSE])))
  ret2 <- ret %>% select(-data)
  multidesign(dsum, ret2, x$column_design)
}

#' @rdname xdata
#' @export
xdata.multidesign <- function(x, ...) x$x

#' @rdname design
#' @export
design.multidesign <- function(x, ...) {
  des <- x$design
  des$.index <- NULL
  des
}

#' @rdname column_design
#' @export
column_design.multidesign <- function(x, ...) {
  x$column_design
}

#' @rdname select_variables
#' @export
select_variables.multidesign <- function(x, ...) {
  cd <- x$column_design
  cd$.col_pos <- seq_len(nrow(cd))
  filt_expr <- rlang::enquos(...)
  cd_filtered <- dplyr::filter(cd, !!!filt_expr)
  if (nrow(cd_filtered) == 0) {
    stop("no columns match the filter expression")
  }
  col_idx <- cd_filtered$.col_pos
  cd_filtered$.col_pos <- NULL
  multidesign(x$x[, col_idx, drop = FALSE], x$design, cd_filtered)
}

#' @rdname fold_over
#' @param preserve_row_ids Logical; if `TRUE`, carry original source row ids into
#'   fold `analysis` and `assessment` designs via a reserved `.orig_index` column.
#'   Matching `held_out$row_ids` metadata is also included.
#' @export
#' @importFrom deflist deflist
fold_over.multidesign <- function(x, ..., preserve_row_ids = FALSE) {
  args <- rlang::enquos(...)
  splits <- split_indices(x, !!!args)
  foldframe <- splits %>% mutate(.fold = seq_len(n()))
  build_multidesign_foldlist(x, foldframe, preserve_row_ids = preserve_row_ids)
}


#' Print Method for Multidesign Objects
#'
#' @description
#' Displays a detailed summary of a multidesign object, including data dimensions,
#' design variables, and column metadata if present.
#'
#' @param x A multidesign object
#' @param ... Additional arguments passed to print methods
#' @return Invisibly returns the input object
#' @method print multidesign
#' @export
print.multidesign <- function(x, ...) {
  # Header
  cat(crayon::bold(crayon::blue("\n=== Multidesign Object ===\n")))

  # Data matrix dimensions
  cat("\nData Matrix: \n")
  cat("  ", nrow(x$x), "observations x", ncol(x$x), "variables \n")

  # Design variables
  cat("\nDesign Variables: \n")
  design_vars <- names(x$design)[!names(x$design) %in% c(".index", ".orig_index")]
  for (var in design_vars) {
    unique_vals <- unique(x$design[[var]])
    n_unique <- length(unique_vals)
    sample_vals <- if (n_unique > 5) {
      paste0(paste(head(unique_vals, 3), collapse=", "),
             "...",
             paste(tail(unique_vals, 2), collapse=", "))
    } else {
      paste(unique_vals, collapse=", ")
    }
    cat("  ", crayon::white("*"), " ", var, ": ",
        crayon::green(n_unique), " levels (", sample_vals, ")\n", sep="")
  }

  # Column metadata if present
  if (!is.null(x$column_design) && ncol(x$column_design) > 0) {
    cat("\nColumn Metadata:\n")
    col_vars <- names(x$column_design)
    for (var in col_vars) {
      unique_vals <- unique(x$column_design[[var]])
      n_unique <- length(unique_vals)
      sample_vals <- if (n_unique > 5) {
        paste0(paste(head(unique_vals, 3), collapse=", "),
               "...",
               paste(tail(unique_vals, 2), collapse=", "))
      } else {
        paste(unique_vals, collapse=", ")
      }
      cat("  ", crayon::white("*"), " ", var, ": ",
          crayon::green(n_unique), " levels (", sample_vals, ")\n", sep="")
    }
  }

  # Footer
  cat(crayon::bold(crayon::blue("\n=======================\n")))
  invisible(x)
}

#' Print Method for Reduced Multidesign Objects
#'
#' @description
#' Displays a detailed summary of a reduced multidesign object, including original
#' and reduced dimensions, design variables, and column metadata if present.
#'
#' @param x A reduced_multidesign object
#' @param ... Additional arguments passed to print methods
#' @return Invisibly returns the input object
#' @method print reduced_multidesign
#' @export
print.reduced_multidesign <- function(x, ...) {
  cat(crayon::bold(crayon::blue("\n=== Reduced Multidesign Object ===\n")))

  # Dimension information
  cat(crayon::bold("\nDimensionality:"), "\n")
  cat("  Original: ", crayon::green(multivarious::shape(x$projector)[1]), " variables\n", sep="")
  cat("  Reduced:  ", crayon::green(multivarious::shape(x$projector)[2]), " components\n", sep="")

  # Data dimensions
  cat(crayon::bold("\nData Matrix:"), "\n")
  cat("  ", crayon::green(paste0(nrow(x$x), " observations x ", ncol(x$x), " components")), "\n")

  # Design variables
  cat(crayon::bold("\nDesign Variables:"), "\n")
  design_vars <- names(x$design)
  design_vars <- design_vars[!design_vars %in% c(".index", ".orig_index")]
  for (var in design_vars) {
    unique_vals <- unique(x$design[[var]])
    n_unique <- length(unique_vals)
    sample_vals <- if (n_unique > 5) {
      paste0(paste(head(unique_vals, 3), collapse=", "),
             "...",
             paste(tail(unique_vals, 2), collapse=", "))
    } else {
      paste(unique_vals, collapse=", ")
    }
    cat("  ", crayon::white("*"), " ", var, ": ",
        crayon::green(n_unique), " levels (", sample_vals, ")\n", sep="")
  }

  # Column metadata if present
  if (!is.null(x$column_design) && ncol(x$column_design) > 0) {
    cat(crayon::bold("\nColumn Metadata:"), "\n")
    col_vars <- names(x$column_design)
    for (var in col_vars) {
      unique_vals <- unique(x$column_design[[var]])
      n_unique <- length(unique_vals)
      sample_vals <- if (n_unique > 5) {
        paste0(paste(head(unique_vals, 3), collapse=", "),
               "...",
               paste(tail(unique_vals, 2), collapse=", "))
      } else {
        paste(unique_vals, collapse=", ")
      }
      cat("  ", crayon::white("*"), " ", var, ": ",
          crayon::green(n_unique), " levels (", sample_vals, ")\n", sep="")
    }
  }

  cat(crayon::bold(crayon::blue("\n=======================\n")))
  invisible(x)
}
