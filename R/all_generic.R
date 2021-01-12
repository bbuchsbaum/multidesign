#' create an multidesign object
#'
#' construct a new multivariate design object linking vector-valued observations and arbitrary design variables
#'
#' @param x the multivariate data (a matrix, a list, or other data container)
#' @param y a design matrix with same number of rows/elements as x
#' @param ... extra args
#' @export
#' @examples
#'
#' X <- matrix(rnorm(20*100), 20, 100)
#' Y <- tibble(condition=rep(letters[1:5], 4))
#'
#' mds <- multidesign(X,Y)
#' sdes <- split(mds, condition)
#' @export
multidesign <- function(x, y, ...) UseMethod("multidesign")

#' create a hyperdesign
#'
#' construct a new `hyperdesign` object. A `hyperdesign` encapusulates a collection of `multidesign` objeects.
#' It is used to model multigroup, muliblock, or multiview datasets, here each block/group/view is
#' associated with a matrix-variate response and an arbitrary design.
#'
#' @param x
#' @param ... extra args
#' @export
hyperdesign <- function(x, ...) UseMethod("hyperdesign")



#' create an multiframe object
#'
#' construct a new multivariate design object linking vector-valued observation and arbitrary design variables
#'
#' @param x the multivariate data (a matrix, a list, or other data container)
#' @param y a design matrix with same number of rows/elements as x
#' @param ... extra args
#' @export
#' @examples
#'
#' X <- matrix(rnorm(20*100), 20, 100)
#' Y <- tibble(condition=rep(letters[1:5], 4))
#'
#' mds <- multidesign(X,Y)
#' sdes <- split(mds, condition)
#' @export
multiframe <- function(x, y, ...) UseMethod("multiframe")


#' summarize data over grouping variable(s)
#'
#' apply a mutlivariate columnwise summary function (e.g. `colMeans`) to subsets
#' of observations formed by one or more grouping variables.
#'
#' @param x the multivariate data (a matrix, a list, or other data container)
#' @param sfun the columnwise summary function (e.g. `colMeans`)
#' @param ... the grouping variables
#'
#'
#' @export
summarize_by <- function(x, sfun, ...) UseMethod("summarize_by")


#' generate cross-validation folds
#'
#' create a set of cross-validation objects over a blocking variable (`stratum`)
#'
#' @param x the dataset to fold over
#' @param ... extra args
#' @export
fold_over <- function(x, ...) UseMethod("fold_over")

#' extract data
#'
#' get the `x` data block
#'
#' @export
xdata <- function(x, ...) UseMethod("xdata")

#' design
#'
#' get design variables
#'
#' @export
design <- function(x, ...) UseMethod("design")


#' split indices
#'
#' extract the row indices of a table, split by one or more variables
#'
#' @param x the object to split
#' @param ... the splitting variables
#' @export
split_indices <- function(x, ...) UseMethod("split_indices")


