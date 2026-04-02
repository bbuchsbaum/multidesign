# Create Cross-Validation Folds from a Hyperdesign Object

Creates cross-validation folds by splitting the data based on specified
design variables. Each fold consists of a training set (analysis) and a
test set (assessment).

## Usage

``` r
# S3 method for class 'hyperdesign'
fold_over(
  x,
  ...,
  inclusion_condition = list(),
  exclusion_condition = list(),
  preserve_row_ids = FALSE
)
```

## Arguments

- x:

  A hyperdesign object

- ...:

  Unquoted names of variables to split on (e.g., condition, subject,
  run)

- inclusion_condition:

  Optional list specifying values to include in assessment sets

- exclusion_condition:

  Optional list specifying values to exclude from assessment sets

- preserve_row_ids:

  Logical; if \`TRUE\`, carry original source row ids into fold
  \`analysis\` and \`assessment\` designs via a reserved \`.orig_index\`
  column. Matching \`held_out\$row_ids\` metadata is also included.

## Value

A foldlist object containing the cross-validation folds

## Details

The function creates folds by splitting the data based on unique
combinations of the specified variables. For each fold, one combination
is held out as the assessment set, while the rest form the analysis set.

Important considerations: \* If a splitting variable is confounded with
blocks (e.g., each subject is in a separate block), the function will
fail as there would be no training data available for that block. \*
Numeric variables (like run numbers) are handled by converting them to
factors for splitting. \* The function preserves the design structure
within each fold.

## Examples

``` r
d1 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(condition = rep(c("A","B"), 5), run = rep(1:5, 2)))
d2 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(condition = rep(c("A","B"), 5), run = rep(1:5, 2)))
hd <- hyperdesign(list(d1, d2))

# Leave-one-block-out folds
folds_block <- fold_over(hd)

# Fold by condition within blocks
folds_cond <- fold_over(hd, condition)
```
