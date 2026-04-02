# Reduce Dimensionality of a Multidesign Object

Performs dimensionality reduction on the data matrix of a multidesign
object while preserving the design structure. By default uses PCA, but
supports any reduction method that returns a projector object.

## Usage

``` r
reduce.multidesign(
  x,
  nc = 2,
  ...,
  rfun = function(x) multivarious::pca(x$x, ncomp = nc, ...)
)
```

## Arguments

- x:

  A multidesign object

- nc:

  Number of components to retain in the reduction

- ...:

  Additional arguments passed to rfun

- rfun:

  Function to perform dimensionality reduction, must return a projector
  object

## Value

A reduced_multidesign object containing:

- x:

  The reduced data matrix

- design:

  Original design information

- projector:

  The projection object used for reduction

## See also

[`multidesign`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md)

Other multidesign functions:
[`multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md),
[`split.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md),
[`split_indices.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split_indices.multidesign.md),
[`subset.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multidesign.md),
[`summarize_by.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multidesign.md)

## Examples

``` r
if (FALSE) { # \dontrun{
X <- matrix(rnorm(100*20), 100, 20)
Y <- tibble::tibble(condition = rep(letters[1:4], each=25))
mds <- multidesign(X, Y)

# Reduce to 5 components using PCA
reduced_mds <- reduce(mds, nc=5)
} # }
```
