# Query an Explicit Cell-Observation Mask

A cell mask records which entries of a multidesign data matrix are
observed. It is independent of the values stored in \`x\`: missing
values in \`x\` never create a mask, and an all-\`FALSE\` mask row
remains a design row.

## Usage

``` r
cell_mask(x, ...)

has_cell_mask(x, ...)

# S3 method for class 'multidesign'
cell_mask(x, ...)

# S3 method for class 'multidesign'
has_cell_mask(x, ...)
```

## Arguments

- x:

  A multidesign object.

- ...:

  Additional arguments passed to methods.

## Value

\`cell_mask()\` returns a logical matrix or \`NULL\` when no explicit
mask is stored. \`has_cell_mask()\` returns a length-one logical value.

## Examples

``` r
md <- multidesign(
  matrix(c(NA, 2, 3, 4), nrow = 2),
  data.frame(group = c("a", "b")),
  cells = matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2)
)
has_cell_mask(md)
#> [1] TRUE
cell_mask(md)
#>       [,1] [,2]
#> [1,]  TRUE TRUE
#> [2,] FALSE TRUE
# The observed NA remains an ordinary stored value.
cell_mask(md)[1, 1]
#> [1] TRUE
```
