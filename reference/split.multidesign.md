# Split a Multidesign Object by Design Variables

Splits a multidesign object into multiple objects based on combinations
of one or more design variables.

## Usage

``` r
# S3 method for class 'multidesign'
split(x, f, drop = FALSE, ...)
```

## Arguments

- x:

  A multidesign object

- f:

  Unquoted name of the first variable to split by

- drop:

  Logical; unused, present for compatibility with the generic

- ...:

  Additional unquoted names of variables to split by

## Value

A list of multidesign objects, one for each combination of splitting
variables

## See also

[`subset.multidesign`](https://bbuchsbaum.github.io/multidesign/reference/subset.multidesign.md)

Other multidesign functions:
[`multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md),
[`reduce.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/reduce.multidesign.md),
[`split_indices.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split_indices.multidesign.md),
[`subset.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multidesign.md),
[`summarize_by.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multidesign.md)

## Examples

``` r
X <- matrix(rnorm(100*20), 100, 20)
Y <- tibble::tibble(
  condition = rep(c("A", "B"), each=50),
  block = rep(1:2, times=50)
)
mds <- multidesign(X, Y)

# Split by condition
split_by_cond <- split(mds, condition)

# Split by condition and block
split_by_both <- split(mds, condition, block)
```
