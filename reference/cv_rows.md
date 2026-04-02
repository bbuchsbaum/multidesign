# Create Cross-validation Folds from Explicit Row Indices

Creates cross-validation folds from explicit assessment row indices
rather than deriving folds from design variables. This is useful when
row holdouts have already been chosen and you want to reuse the
package's existing \`foldlist\` and \`cross_validate()\` machinery.

## Usage

``` r
cv_rows(x, rows, ...)

# S3 method for class 'multidesign'
cv_rows(x, rows, preserve_row_ids = FALSE, ...)

# S3 method for class 'multiframe'
cv_rows(x, rows, preserve_row_ids = FALSE, ...)

# S3 method for class 'hyperdesign'
cv_rows(x, rows, preserve_row_ids = FALSE, ...)
```

## Arguments

- x:

  The dataset to fold over (\`multidesign\`, \`hyperdesign\`, or
  \`multiframe\`).

- rows:

  Assessment row indices for each fold. \* For \`multidesign\` and
  \`multiframe\`, supply a list of integer vectors, one per fold. \* For
  \`hyperdesign\`, supply a list of per-fold block mappings. Each fold
  is a list whose elements are integer vectors of row indices, keyed by
  block name or block position.

- ...:

  Additional arguments passed to methods (currently unused).

- preserve_row_ids:

  Logical; if \`TRUE\`, carry original source row ids into fold
  \`analysis\` and \`assessment\` designs via a reserved \`.orig_index\`
  column. Where \`held_out\` metadata is present, matching \`row_ids\`
  are also included.

## Value

A \`foldlist\` object containing \`analysis\`, \`assessment\`, and
optional \`held_out\` metadata.

## See also

[`fold_over`](https://bbuchsbaum.github.io/multidesign/reference/fold_over.md)
for creating folds from design variables,
[`cross_validate`](https://bbuchsbaum.github.io/multidesign/reference/cross_validate.md)
for executing the fit/score loop

## Examples

``` r
X <- matrix(rnorm(40), 10, 4)
Y <- data.frame(group = rep(c("A", "B"), each = 5))
mds <- multidesign(X, Y)

folds <- cv_rows(mds, rows = list(1:2, 6:7))
folds[[1]]$assessment
#> 
#> === Multidesign Object ===
#> 
#> Data Matrix: 
#>    2 observations x 4 variables 
#> 
#> Design Variables: 
#>   * group: 1 levels (A)
#> 
#> Column Metadata:
#>   * .index: 4 levels (1, 2, 3, 4)
#> 
#> =======================
#> 
folds_with_ids <- cv_rows(mds, rows = list(1:2), preserve_row_ids = TRUE)
folds_with_ids[[1]]$assessment$design$.orig_index
#> [1] 1 2

d1 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(run = 1:6))
d2 <- multidesign(matrix(rnorm(30), 6, 5), data.frame(run = 1:6))
hd <- hyperdesign(list(d1, d2), block_names = c("subj1", "subj2"))

hd_folds <- cv_rows(hd, rows = list(
  list(subj1 = 1:2, subj2 = 1:2),
  list(subj1 = 3:4, subj2 = 3:4)
))
hd_folds[[1]]$assessment
#> 
#> === Hyperdesign Object ===
#> 
#> Number of blocks:  2 
#> 
#> +- Block  1  (subj1)  -----------------
#> |  Dimensions: 2 x 5 
#> |  Design Variables: run 
#> |  Design Structure: 
#> |   * run: 2 levels (1, 2)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> +- Block  2  (subj2)  -----------------
#> |  Dimensions: 2 x 5 
#> |  Design Variables: run 
#> |  Design Structure: 
#> |   * run: 2 levels (1, 2)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> =======================
#> 
```
