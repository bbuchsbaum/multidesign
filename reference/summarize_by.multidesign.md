# Summarize a Multidesign Object by Design Variables

Computes summaries of the data matrix grouped by combinations of design
variables.

## Usage

``` r
# S3 method for class 'multidesign'
summarize_by(x, ..., sfun = colMeans, extract_data = FALSE, aggregate = NULL)
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

- aggregate:

  Optional explicit cellwise aggregation rule. Required when \`x\` has a
  cell mask; accepts \`"mean"\` or a scalar-returning function.

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

# Masked summaries require an explicit coordinate-wise rule
mask <- matrix(TRUE, nrow(X), ncol(X))
mask[1:10, 1] <- FALSE
masked <- multidesign(X, Y, cells = mask)
masked_means <- summarize_by(masked, condition, aggregate = "mean")
cell_mask(masked_means)
#>         [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> group 1 TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE  TRUE  TRUE  TRUE  TRUE
#> group 2 TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE  TRUE  TRUE  TRUE  TRUE
#>         [,14] [,15] [,16] [,17] [,18] [,19] [,20]
#> group 1  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> group 2  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
```
