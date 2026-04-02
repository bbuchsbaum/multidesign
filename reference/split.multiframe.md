# Split a Multiframe by Design Variables

Splits a multiframe into groups based on one or more design variables.
Optionally collapses the observations within each group into a single
matrix.

## Usage

``` r
# S3 method for class 'multiframe'
split(x, ..., collapse = FALSE)
```

## Arguments

- x:

  A multiframe object

- ...:

  Unquoted names of variables to split by

- collapse:

  Logical; if TRUE, combines observations in each group into a matrix

## Value

A nested tibble containing the split data

## See also

Other multiframe functions:
[`[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-.observation_set.md),
[`[[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-sub-.observation_set.md),
[`multiframe.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.list.md),
[`multiframe.matrix()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.matrix.md),
[`obs_group()`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md),
[`subset.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multiframe.md),
[`summarize_by.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multiframe.md),
[`xdata.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.multiframe.md)

## Examples

``` r
x <- matrix(1:40, 10, 4)
y <- data.frame(
  condition = rep(c("A", "B"), each=5),
  subject = rep(1:5, times=2)
)
mf <- multiframe(x, y)

# Split by condition
split_by_cond <- split(mf, condition)

# Split by condition and subject, collapsing observations
split_by_both <- split(mf, condition, subject, collapse=TRUE)
```
