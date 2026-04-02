# Split Indices by Variables

Extracts the row indices of a dataset, split by one or more variables.
This function is useful for creating grouped analyses, stratified
sampling, or preparing data for cross-validation without actually
duplicating the data.

## Usage

``` r
split_indices(x, ...)

# S3 method for class 'multiframe'
split_indices(x, ..., collapse = FALSE)
```

## Arguments

- x:

  The object to split (data frame, multidesign, multiframe, etc.)

- ...:

  Unquoted names of variables to split by

- collapse:

  Logical; whether to collapse the resulting indices

## Value

A tibble or list containing:

- group variables:

  The splitting variables and their values

- indices:

  List column of integer vectors with row indices for each group

## Details

The function returns indices rather than the actual data, which is
memory-efficient for large datasets. It works with various object types:
\* For data frames: Returns indices grouped by specified variables \*
For multidesign objects: Returns indices grouped by design variables \*
For multiframe objects: Returns indices grouped by design variables

The resulting indices can be used to subset the original data or to
create training/testing splits for machine learning.

## See also

[`fold_over`](https://bbuchsbaum.github.io/multidesign/reference/fold_over.md)
for creating cross-validation folds,
[`split`](https://rdrr.io/r/base/split.html) for splitting the actual
data

## Examples

``` r
# With a multidesign object
mds <- multidesign(matrix(rnorm(100*10), 100, 10),
                  data.frame(group=rep(letters[1:4], each=25),
                             block=rep(1:5, times=20)))
mds_indices <- split_indices(mds, group)
```
