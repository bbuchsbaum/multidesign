# Generate Cross-validation Folds

Creates cross-validation folds by splitting data based on specified
design variables or stratification factors. This function helps in
creating stratified cross-validation splits while respecting the
structure of the data, particularly for complex experimental designs.

## Usage

``` r
fold_over(x, ...)

# S3 method for class 'multidesign'
fold_over(x, ..., preserve_row_ids = FALSE)

# S3 method for class 'multiframe'
fold_over(x, ..., preserve_row_ids = FALSE)
```

## Arguments

- x:

  The dataset to fold over (multidesign, hyperdesign, matrix, etc.)

- ...:

  Additional arguments passed to methods: \* For
  multidesign/hyperdesign: unquoted names of variables to split on \*
  For multidesign/hyperdesign/multiframe: \`preserve_row_ids = TRUE\` to
  carry original source row ids into fold outputs via a reserved
  \`.orig_index\` column \* For matrices: nfolds (number of folds),
  stratum (stratification variable)

- preserve_row_ids:

  Logical; if \`TRUE\`, carry original source row ids into fold
  \`analysis\` and \`assessment\` designs via a reserved \`.orig_index\`
  column. Matching \`held_out\$row_ids\` metadata is also included.

## Value

A foldlist object containing:

- analysis:

  Training data for each fold

- assessment:

  Test data for each fold

- held_out:

  Information about which values were held out in each fold

## Details

The function creates folds by splitting the data based on unique
combinations of the specified variables. For each fold, one combination
is held out as the assessment set (test set), while the rest form the
analysis set (training set).

Different behaviors are implemented for different object types: \* For
multidesign objects: Creates folds based on design variables \* For
hyperdesign objects: Creates folds across blocks, respecting block
structure \* For matrices: Creates folds based on external
stratification variables

## See also

[`multidesign`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md)
for creating multidesign objects,
[`hyperdesign`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.md)
for creating hyperdesign objects,
[`split_indices`](https://bbuchsbaum.github.io/multidesign/reference/split_indices.md)
for splitting indices without creating folds

## Examples

``` r
# With a multidesign object
X <- matrix(rnorm(100*10), 100, 10)
groups <- rep(1:2, each=50)
mds <- multidesign(X, data.frame(group=groups, subject=rep(1:10, each=10)))
folds_by_group <- fold_over(mds, group)

# With a hyperdesign object (multiple subjects)
d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
                 data.frame(condition=rep(c("A","B"), 5), run=1:10))
d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
                 data.frame(condition=rep(c("A","B"), 5), run=1:10))
hd <- hyperdesign(list(d1, d2), block_names=c("subject1", "subject2"))
folds_by_condition <- fold_over(hd, condition)
```
