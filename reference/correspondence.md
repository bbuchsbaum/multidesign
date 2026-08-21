# Query Row Correspondence in a Hyperdesign

A correspondence object describes how local rows in each hyperdesign
block map to a global entity universe. Correspondence is opt-in:
construct the hyperdesign with either \`id\` or \`positional = TRUE\`.

## Usage

``` r
correspondence(x, ...)

# S3 method for class 'hyperdesign'
correspondence(x, sort_ids = FALSE, ...)
```

## Arguments

- x:

  A hyperdesign object.

- ...:

  Additional arguments passed to methods.

- sort_ids:

  Logical; if \`TRUE\`, order normalized global identifiers with base R
  radix ordering. The default preserves first-seen order.

## Value

For a contracted hyperdesign, a \`md_correspondence\` list containing
the declared key, assumption, global identifiers, local identifiers,
integer row maps, incidence matrix, pairwise overlap counts, and block
connectivity flag.

## Details

Explicit identifiers are compared through a common character
representation, so integer \`1\` and character \`"1"\` match.
\`global_ids\` follows first-seen block and row order, and each integer
\`row_map\` indexes that vector. Missing values in the data matrix are
ordinary data and are never interpreted as a correspondence mask.

## Examples

``` r
a <- multidesign(matrix(1:6, ncol = 2), data.frame(id = c("B", "A", "C")))
b <- multidesign(matrix(7:10, ncol = 2), data.frame(id = c("A", "D")))
hd <- hyperdesign(list(a = a, b = b), id = "id", space = "common")
correspondence(hd)$global_ids
#> [1] "B" "A" "C" "D"
correspondence(hd, sort_ids = TRUE)$global_ids
#> [1] "A" "B" "C" "D"
```
