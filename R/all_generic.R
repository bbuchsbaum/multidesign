#' @importFrom crayon bold blue green white yellow
NULL

#' create an observation object
#'
#' construct a new multivariate observation vector
#'
#' @param x the data source
#' @param i the index of the observation
#' @export
observation <- function(x, i) UseMethod("observation")

# ... rest of existing file content remains unchanged ...
