# Summarize a Multiframe by Design Variables

Computes summaries of observations grouped by design variables.

## Usage

``` r
# S3 method for class 'multiframe'
summarize_by(x, ..., sfun = colMeans, extract_data = FALSE)
```

## Arguments

- x:

  A multiframe object

- ...:

  Unquoted names of variables to group by

- sfun:

  Summary function to apply (default is colMeans)

- extract_data:

  Logical; if TRUE, returns just the summarized data

## Value

A tibble containing summarized data for each group

## See also

Other multiframe functions:
[`[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-.observation_set.md),
[`[[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-sub-.observation_set.md),
[`multiframe.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.list.md),
[`multiframe.matrix()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.matrix.md),
[`obs_group()`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md),
[`split.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/split.multiframe.md),
[`subset.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multiframe.md),
[`xdata.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.multiframe.md)

## Examples

``` r
x <- matrix(1:40, 10, 4)
y <- data.frame(
  condition = rep(c("A", "B"), each=5),
  subject = rep(1:5, times=2)
)
mf <- multiframe(x, y)

# Get means by condition
means_by_cond <- summarize_by(mf, condition)
```
