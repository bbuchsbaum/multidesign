# Get Split Indices for a Multidesign Object

Returns indices for splitting a multidesign object based on combinations
of design variables, without actually splitting the data.

## Usage

``` r
# S3 method for class 'multidesign'
split_indices(x, ..., collapse = FALSE)
```

## Arguments

- x:

  A multidesign object

- ...:

  Unquoted names of variables to split by

- collapse:

  Logical; whether to collapse the resulting indices

## Value

A tibble containing:

- group variables:

  The splitting variables and their values

- indices:

  List column of indices for each group

- .splitvar:

  Combined string of grouping values

## See also

[`split.multidesign`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md)

Other multidesign functions:
[`multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md),
[`reduce.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/reduce.multidesign.md),
[`split.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md),
[`subset.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multidesign.md),
[`summarize_by.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multidesign.md)

## Examples

``` r
X <- matrix(rnorm(40*10), 40, 10)
Y <- data.frame(condition = rep(c("A","B"), each=20),
                block = rep(1:4, each=10))
mds <- multidesign(X, Y)
idx <- split_indices(mds, condition)
```
