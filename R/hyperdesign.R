#' Subset a Hyperdesign Object
#'
#' Create a new hyperdesign object containing only observations that meet specified
#' conditions based on design variables.
#'
#' @param x A hyperdesign object
#' @param subset Expression used to filter observations based on design variables
#' @param ... Additional arguments (not used)
#' @return A new hyperdesign object containing only the selected observations
#' @family hyperdesign functions
#' @method subset hyperdesign
#' @export
#' @examples
#' # Create example hyperdesign
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(subject=1, condition=rep(c("A","B"), 5)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(subject=2, condition=rep(c("A","B"), 5)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(subject=3, condition=rep(c("A","B"), 5)))
#' hd <- hyperdesign(list(d1, d2, d3))
#'
#' # Keep only observations where condition == "A"
#' subset_hd <- subset(hd, condition == "A")
subset.hyperdesign <- function(x, subset, ...) {
  out <- lapply(x, function(d) {
    subset(d, !!rlang::enquo(subset))
  })

  rem <- unlist(purrr::map(out, is.null))
  if (sum(rem) == length(x)) {
    stop("subset expression does not match any rows in hyperdesign `x`")
  }

  onam <- names(x)[!rem]
  hyperdesign(out[!rem], onam)
}
