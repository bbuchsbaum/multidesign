# Summarize a Multidesign Object by Design Variables

Computes summaries of the data matrix grouped by combinations of design
variables.

## Usage

``` r
# S3 method for class 'multidesign'
summarize_by(x, ..., sfun = colMeans, extract_data = FALSE)
```

## Arguments

- x:

  A multidesign object

- ...:

  Unquoted names of variables to group by

- sfun:

  Summary function to apply (default is colMeans)

- extract_data:

  Logical; whether to extract raw data instead of computing summary

## Value

A new multidesign object containing:

- x:

  Matrix of summary statistics

- design:

  Design information for each summary

- column_design:

  Original column metadata

## See also

[`split.multidesign`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md)

Other multidesign functions:
[`multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md),
[`reduce.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/reduce.multidesign.md),
[`split.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md),
[`split_indices.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split_indices.multidesign.md),
[`subset.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multidesign.md)

## Examples

``` r
X <- matrix(rnorm(100*20), 100, 20)
Y <- tibble::tibble(
  condition = rep(c("A", "B"), each=50),
  block = rep(1:2, times=50)
)
mds <- multidesign(X, Y)

# Get means by condition
means_by_cond <- summarize_by(mds, condition)

# Get means by condition and block
means_by_both <- summarize_by(mds, condition, block)
```
