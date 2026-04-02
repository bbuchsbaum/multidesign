# Subset a Multiframe Object

Creates a new multiframe object containing only observations that meet
specified conditions based on design variables. The lazy-evaluated
observation closures remain valid in the subset.

## Usage

``` r
# S3 method for class 'multiframe'
subset(x, fexpr, ...)
```

## Arguments

- x:

  A multiframe object

- fexpr:

  An expression used to filter observations based on design variables

- ...:

  Additional arguments (not used)

## Value

A new multiframe object containing only matching observations, or NULL
if no matches

## See also

Other multiframe functions:
[`[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-.observation_set.md),
[`[[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-sub-.observation_set.md),
[`multiframe.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.list.md),
[`multiframe.matrix()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.matrix.md),
[`obs_group()`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md),
[`split.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/split.multiframe.md),
[`summarize_by.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multiframe.md),
[`xdata.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.multiframe.md)

## Examples

``` r
X <- matrix(rnorm(20 * 5), 20, 5)
Y <- data.frame(
  condition = rep(c("A", "B"), each = 10),
  block = rep(1:4, each = 5)
)
mf <- multiframe(X, Y)
mf_A <- subset(mf, condition == "A")
```
