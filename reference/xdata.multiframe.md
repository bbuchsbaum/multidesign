# Extract Data Matrix from a Multiframe Object

Materializes all lazy-evaluated observations into a single data matrix.
Each observation is evaluated and the results are row-bound together.

## Usage

``` r
# S3 method for class 'multiframe'
xdata(x, ...)
```

## Arguments

- x:

  A multiframe object

- ...:

  Additional arguments (not used)

## Value

A matrix containing all materialized observations

## See also

Other multiframe functions:
[`[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-.observation_set.md),
[`[[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-sub-.observation_set.md),
[`multiframe.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.list.md),
[`multiframe.matrix()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.matrix.md),
[`obs_group()`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md),
[`split.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/split.multiframe.md),
[`subset.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multiframe.md),
[`summarize_by.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multiframe.md)

## Examples

``` r
X <- matrix(rnorm(20 * 5), 20, 5)
Y <- data.frame(condition = rep(c("A", "B"), each = 10))
mf <- multiframe(X, Y)
mat <- xdata(mf)  # Returns the full 20x5 matrix
```
