# Subset a Multidesign Object

Creates a new multidesign object containing only observations that meet
specified conditions based on design variables.

## Usage

``` r
# S3 method for class 'multidesign'
subset(x, fexpr, ...)
```

## Arguments

- x:

  A multidesign object

- fexpr:

  An expression used to filter observations based on design variables

- ...:

  Additional arguments (not used)

## Value

A new multidesign object containing only matching observations, or NULL
if no matches

## See also

[`split.multidesign`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md)

Other multidesign functions:
[`multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md),
[`reduce.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/reduce.multidesign.md),
[`split.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md),
[`split_indices.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split_indices.multidesign.md),
[`summarize_by.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multidesign.md)

## Examples

``` r
X <- matrix(rnorm(100*20), 100, 20)
Y <- tibble::tibble(
  condition = rep(c("A", "B"), each=50),
  block = rep(1:2, times=50)
)
mds <- multidesign(X, Y)

# Subset to condition A only
mds_A <- subset(mds, condition == "A")

# Subset to block 1 and condition A
mds_A1 <- subset(mds, condition == "A" & block == 1)
```
