#' Execute Cross-Validation Over Folds
#'
#' @description
#' Runs a cross-validation pipeline by iterating over folds, fitting a model
#' on the analysis (training) set, and scoring predictions on the assessment
#' (test) set. Returns a structured result object for summarizing performance.
#'
#' @param folds A foldlist object (e.g., from \code{\link{fold_over}})
#' @param fit_fn A function that takes the analysis set and returns a model object.
#'   Signature: \code{fit_fn(analysis)} where analysis is a multidesign, hyperdesign,
#'   or multiframe object.
#' @param score_fn A function that takes a model and the assessment set, and returns
#'   a single numeric value, a named numeric vector, or a named list of numeric values.
#'   Signature: \code{score_fn(model, assessment)}
#' @param ... Additional arguments (currently unused)
#'
#' @return A \code{cv_result} object containing a tibble of per-fold scores
#'
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' Y <- data.frame(group = rep(1:5, each = 20))
#' mds <- multidesign(X, Y)
#' folds <- fold_over(mds, group)
#'
#' result <- cross_validate(
#'   folds,
#'   fit_fn = function(analysis) {
#'     list(mean = colMeans(xdata(analysis)))
#'   },
#'   score_fn = function(model, assessment) {
#'     pred_error <- mean((xdata(assessment) - matrix(model$mean, nrow = nrow(xdata(assessment)),
#'                         ncol = length(model$mean), byrow = TRUE))^2)
#'     c(mse = pred_error)
#'   }
#' )
#' print(result)
#' summary(result)
#'
#' @seealso
#'   \code{\link{fold_over}} for creating folds,
#'   \code{\link{new_cv_result}} for the result structure
#' @export
cross_validate <- function(folds, fit_fn, score_fn, ...) {
  chk::chk_s3_class(folds, "foldlist")
  chk::chk_function(fit_fn)
  chk::chk_function(score_fn)

  all_scores <- lapply(seq_along(folds), function(i) {
    fold <- folds[[i]]
    model <- fit_fn(fold$analysis)
    score <- score_fn(model, fold$assessment)

    # Normalize score to a named list
    if (is.list(score)) {
      if (!all(sapply(score, is.numeric))) {
        stop("score_fn must return numeric values")
      }
      score_df <- as.data.frame(score, stringsAsFactors = FALSE)
    } else if (is.numeric(score)) {
      if (is.null(names(score))) {
        score_df <- data.frame(score = score)
      } else {
        score_df <- as.data.frame(as.list(score), stringsAsFactors = FALSE)
      }
    } else {
      stop("score_fn must return a numeric value, named numeric vector, or named list")
    }

    score_df$.fold <- i
    score_df
  })

  scores <- do.call(rbind, all_scores)
  new_cv_result(scores, call = match.call())
}
