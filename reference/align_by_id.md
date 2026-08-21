# Align Hyperdesign Blocks by Their Entity Contract

Places every common-space block on the global entity universe described
by \[correspondence()\]. The result keeps row presence and cell
observation as separate logical arrays; the fill value is presentation
only.

## Usage

``` r
align_by_id(x, ...)

# S3 method for class 'hyperdesign'
align_by_id(x, fill = NA_real_, sort_ids = FALSE, ...)
```

## Arguments

- x:

  A hyperdesign with declared correspondence and \`space = "common"\`.

- ...:

  Additional arguments passed to methods.

- fill:

  Length-one numeric value used to display absent rows and explicitly
  unobserved cells.

- sort_ids:

  Logical; if \`TRUE\`, use sorted global identifiers rather than
  first-seen order.

## Value

An \`aligned_hyperdesign\` list containing \`x\`, \`cells\`,
\`observed\`, \`correspondence\`, \`column_design\`, and
\`block_names\`.

## Examples

``` r
a <- multidesign(
  matrix(c(NA, 0, 1, 0), ncol = 2, byrow = TRUE),
  data.frame(id = c("A", "B")),
  cells = matrix(c(TRUE, TRUE, TRUE, FALSE), ncol = 2, byrow = TRUE)
)
b <- multidesign(
  matrix(c(0, 1, 1, 1), ncol = 2, byrow = TRUE),
  data.frame(id = c("B", "C"))
)
hd <- hyperdesign(list(a = a, b = b), id = "id", space = "common")
aligned <- align_by_id(hd)
aligned$observed
#>       a     b
#> A  TRUE FALSE
#> B  TRUE  TRUE
#> C FALSE  TRUE
aligned$cells[, , "a"]
#>       variable
#> entity     1     2
#>      A  TRUE  TRUE
#>      B  TRUE FALSE
#>      C FALSE FALSE
aligned$x[, , "a"]
#>       variable
#> entity  1  2
#>      A NA  0
#>      B  1 NA
#>      C NA NA
```
