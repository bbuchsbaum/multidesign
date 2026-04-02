# Create a Multiframe Object

Constructs a new multivariate design object linking vector-valued
observations with arbitrary design variables. A multiframe combines
experimental design metadata with lazy-evaluated observations, providing
a flexible interface for data manipulation.

## Usage

``` r
multiframe(x, y, ...)
```

## Arguments

- x:

  The multivariate data (a matrix, list, or other data container)

- y:

  A design matrix or data frame with same number of rows/elements as x

- ...:

  Additional arguments passed to methods

## Value

A multiframe object containing:

- design:

  A tibble with design variables and observation functions

## Details

A multiframe object is similar to a multidesign but with enhanced
support for data frame operations. Key features include: \* Lazy
evaluation of observations (only loaded when needed) \* Integration with
tidyverse functions for data manipulation \* Support for various data
sources (matrices, lists, vectors) \* Methods for splitting,
summarizing, and transforming data

## See also

[`multidesign`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md)
for an alternative implementation,
[`observation`](https://bbuchsbaum.github.io/multidesign/reference/observation.md)
for the underlying observation structure,
[`obs_group`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md)
for creating groups of observations

## Examples

``` r
# Create sample data
X <- matrix(rnorm(20*100), 20, 100)
Y <- tibble::tibble(condition = rep(letters[1:5], 4))

# Create multiframe object
mf <- multiframe(X, Y)

# Access first observation
obs1 <- mf$design$.obs[[1]]()

# Split by condition
split_by_cond <- split(mf, condition)

# Summarize by condition
means_by_cond <- summarize_by(mf, condition)
```
