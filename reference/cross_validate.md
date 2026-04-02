# Execute Cross-Validation Over Folds

Runs a cross-validation pipeline by iterating over folds, fitting a
model on the analysis (training) set, and scoring predictions on the
assessment (test) set. Returns a structured result object for
summarizing performance.

## Usage

``` r
cross_validate(folds, fit_fn, score_fn, ...)
```

## Arguments

- folds:

  A foldlist object (e.g., from
  [`fold_over`](https://bbuchsbaum.github.io/multidesign/reference/fold_over.md))

- fit_fn:

  A function that takes the analysis set and returns a model object.
  Signature: `fit_fn(analysis)` where analysis is a multidesign,
  hyperdesign, or multiframe object.

- score_fn:

  A function that takes a model and the assessment set, and returns a
  single numeric value, a named numeric vector, or a named list of
  numeric values. Signature: `score_fn(model, assessment)`. If
  `score_fn` declares `fold`, `fold_id`, or `...`, the current fold
  object and fold index are also supplied.

- ...:

  Additional arguments (currently unused)

## Value

A `cv_result` object containing a tibble of per-fold scores

## See also

[`fold_over`](https://bbuchsbaum.github.io/multidesign/reference/fold_over.md)
for creating folds,
[`new_cv_result`](https://bbuchsbaum.github.io/multidesign/reference/new_cv_result.md)
for the result structure

## Examples

``` r
X <- matrix(rnorm(100 * 10), 100, 10)
Y <- data.frame(group = rep(1:5, each = 20))
mds <- multidesign(X, Y)
folds <- fold_over(mds, group)

result <- cross_validate(
  folds,
  fit_fn = function(analysis) {
    list(mean = colMeans(xdata(analysis)))
  },
  score_fn = function(model, assessment) {
    pred_error <- mean((xdata(assessment) - matrix(model$mean, nrow = nrow(xdata(assessment)),
                        ncol = length(model$mean), byrow = TRUE))^2)
    c(mse = pred_error)
  }
)
print(result)
#> 
#> === Cross-Validation Result ===
#> 
#> Folds: 5 
#> 
#> Metrics:
#>    mse : mean = 0.9542 , sd = 0.0929 
#> 
#> ===============================
#> 
summary(result)
#> # A tibble: 1 × 6
#>   metric  mean     sd   min   max median
#>   <chr>  <dbl>  <dbl> <dbl> <dbl>  <dbl>
#> 1 mse    0.954 0.0929 0.801  1.03  0.982

result_with_context <- cross_validate(
  folds,
  fit_fn = function(analysis) list(mean = colMeans(xdata(analysis))),
  score_fn = function(model, assessment, fold_id) {
    c(mse = mean((xdata(assessment) - model$mean)^2), fold_id = fold_id)
  }
)
```
