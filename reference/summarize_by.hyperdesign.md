# Summarize Blocks of a Hyperdesign

Applies \[summarize_by()\] to each block. Explicit entity correspondence
is preserved only when the entity key remains among the grouping
variables. Otherwise the result is intentionally uncontracted with
respect to rows. Masked blocks require an explicit coordinate-wise
\`aggregate\` rule.

## Usage

``` r
# S3 method for class 'hyperdesign'
summarize_by(x, ..., sfun = colMeans, extract_data = FALSE, aggregate = NULL)
```

## Arguments

- x:

  A hyperdesign object.

- ...:

  Unquoted grouping variables.

- sfun:

  Legacy matrix summary function used for unmasked blocks.

- extract_data:

  Logical legacy argument forwarded to block methods.

- aggregate:

  Explicit cellwise aggregation rule required for masked blocks.

## Value

A summarized hyperdesign.

## Examples

``` r
md <- multidesign(
  matrix(1:8, nrow = 4),
  data.frame(entity = c("A", "A", "B", "B")),
  cells = matrix(c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE),
                 nrow = 4)
)
hd <- hyperdesign(list(one = md))
out <- summarize_by(hd, entity, aggregate = "mean")
cell_mask(out[[1]])
#>         [,1]  [,2]
#> group 1 TRUE  TRUE
#> group 2 TRUE FALSE
```
