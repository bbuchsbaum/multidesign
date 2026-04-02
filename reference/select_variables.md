# Select Variables Based on Column Design

Subsets the columns (variables) of an object based on conditions applied
to the column design metadata. This allows filtering variables by their
properties (e.g., region, type, hemisphere) using dplyr-style filter
expressions.

## Usage

``` r
select_variables(x, ...)

# S3 method for class 'hyperdesign'
select_variables(x, ...)

# S3 method for class 'multidesign'
select_variables(x, ...)
```

## Arguments

- x:

  The object to select variables from (multidesign, hyperdesign, etc.)

- ...:

  Filter expressions evaluated against the column design (using
  dplyr::filter semantics)

## Value

An object of the same class with a subset of columns/variables

## See also

[`column_design`](https://bbuchsbaum.github.io/multidesign/reference/column_design.md)
for extracting column metadata,
[`multidesign`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md)
for creating multidesign objects

## Examples

``` r
X <- matrix(rnorm(20*10), 20, 10)
Y <- data.frame(condition = rep(c("A", "B"), each=10))
col_info <- data.frame(
  region = rep(c("frontal", "parietal"), 5),
  hemisphere = rep(c("left", "right"), each=5)
)
mds <- multidesign(X, Y, col_info)

# Select only frontal variables
mds_frontal <- select_variables(mds, region == "frontal")
```
