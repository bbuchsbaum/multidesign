#' @export
xdata.hyperdesign <- function(x, ...) {
  if (missing(block)) {
    lapply(x, xdata)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, c(1, length(x)))
    xdata(x[[block]])
  }
}
