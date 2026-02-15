#' Create a cv_result Object
#'
#' @description
#' Internal constructor for cv_result objects. Validates that scores is a
#' data frame with a `.fold` column and wraps it in a cv_result structure.
#'
#' @param scores A data frame containing cross-validation scores with a `.fold` column
#' @param call Optional call object recording how the result was created
#' @return A cv_result object
#' @keywords internal
#' @export
new_cv_result <- function(scores, call = NULL) {
  chk::chk_s3_class(scores, "data.frame")
  if (!".fold" %in% names(scores)) {
    stop("scores must contain a '.fold' column")
  }
  scores <- tibble::as_tibble(scores)
  structure(list(scores = scores, call = call), class = "cv_result")
}

#' Print Method for cv_result Objects
#'
#' @description
#' Displays fold count and mean/sd for each numeric score column.
#'
#' @param x A cv_result object
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the input object
#' @importFrom stats sd median
#' @method print cv_result
#' @export
print.cv_result <- function(x, ...) {
  scores <- x$scores
  n_folds <- length(unique(scores$.fold))
  cat(crayon::bold(crayon::blue("\n=== Cross-Validation Result ===\n")))
  cat("\nFolds:", crayon::green(n_folds), "\n")

  numeric_cols <- setdiff(names(scores), ".fold")
  numeric_cols <- numeric_cols[sapply(numeric_cols, function(nm) is.numeric(scores[[nm]]))]

  if (length(numeric_cols) > 0) {
    cat(crayon::bold("\nMetrics:\n"))
    for (nm in numeric_cols) {
      vals <- scores[[nm]]
      cat("  ", nm, ": mean =", crayon::green(round(mean(vals), 4)),
          ", sd =", crayon::green(round(stats::sd(vals), 4)), "\n")
    }
  }

  cat(crayon::bold(crayon::blue("\n===============================\n")))
  invisible(x)
}

#' Summary Method for cv_result Objects
#'
#' @description
#' Returns a tibble with metric/mean/sd/min/max/median per numeric score column.
#'
#' @param object A cv_result object
#' @param ... Additional arguments (not used)
#' @return A tibble with summary statistics for each metric
#' @method summary cv_result
#' @export
summary.cv_result <- function(object, ...) {
  scores <- object$scores
  numeric_cols <- setdiff(names(scores), ".fold")
  numeric_cols <- numeric_cols[sapply(numeric_cols, function(nm) is.numeric(scores[[nm]]))]

  rows <- lapply(numeric_cols, function(nm) {
    vals <- scores[[nm]]
    tibble::tibble(
      metric = nm,
      mean = mean(vals),
      sd = stats::sd(vals),
      min = min(vals),
      max = max(vals),
      median = stats::median(vals)
    )
  })
  do.call(rbind, rows)
}
