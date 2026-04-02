# Create a Multidesign Object

Constructs a new multivariate design object linking vector-valued
observations with design variables. A multidesign object maintains the
relationship between experimental data (observations) and metadata about
experimental conditions (design variables).

Creates a multidesign object that combines experimental data (as a
matrix) with design information (as a data frame) and optional column
metadata. This structure is particularly useful for experimental designs
where observations have multiple associated factors and variables may
have metadata.

## Usage

``` r
multidesign(x, y, ...)

# S3 method for class 'matrix'
multidesign(x, y, column_design = NULL, ...)
```

## Arguments

- x:

  A numeric matrix where rows are observations and columns are variables

- y:

  A data frame containing design variables for each observation (must
  have same number of rows as x)

- ...:

  Additional arguments passed to methods, such as column_design

- column_design:

  Optional data frame containing metadata for columns in x (must have
  same number of rows as ncol(x))

## Value

A multidesign object with components:

- x:

  The input data matrix

- design:

  A tibble containing design variables

- column_design:

  A tibble containing column metadata (if provided)

A multidesign object with components:

- x:

  The input data matrix

- design:

  A tibble containing design variables with an added .index column

- column_design:

  A tibble containing column metadata (empty if not provided)

## Details

A multidesign object consists of three main components: \* A data matrix
where rows represent observations and columns represent variables \* A
design data frame containing experimental factors and conditions for
each observation \* Optional column metadata describing properties of
each variable

A multidesign object consists of three main components: \* A data matrix
where rows represent observations and columns represent variables \* A
design data frame containing experimental factors and conditions for
each observation \* Optional column metadata describing properties of
each variable

The object maintains the relationship between these components while
providing methods for manipulation, subsetting, and analysis.

## See also

[`reduce.multidesign`](https://bbuchsbaum.github.io/multidesign/reference/reduce.multidesign.md)
for dimensionality reduction,
[`split.multidesign`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md)
for splitting by design variables,
[`multiframe`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.md)
for an alternative implementation

[`reduce.multidesign`](https://bbuchsbaum.github.io/multidesign/reference/reduce.multidesign.md)
for dimensionality reduction,
[`split.multidesign`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md)
for splitting by design variables

Other multidesign functions:
[`reduce.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/reduce.multidesign.md),
[`split.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split.multidesign.md),
[`split_indices.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/split_indices.multidesign.md),
[`subset.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multidesign.md),
[`summarize_by.multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multidesign.md)

## Examples

``` r
# Create example data matrix
X <- matrix(rnorm(20*100), 20, 100)

# Create design information
Y <- tibble::tibble(condition=rep(letters[1:5], 4))

# Create multidesign object
mds <- multidesign(X, Y)

# Split by condition
sdes <- split(mds, condition)

# Create example data matrix
X <- matrix(rnorm(20*100), 20, 100)

# Create design information
Y <- tibble::tibble(
  condition = rep(c("control", "treatment"), each=10),
  subject = rep(1:5, times=4)
)

# Create column metadata
col_info <- data.frame(
  roi = paste0("region_", 1:100),
  hemisphere = rep(c("left", "right"), 50)
)

# Create multidesign object
mds <- multidesign(X, Y, col_info)
```
