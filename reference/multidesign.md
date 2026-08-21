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
multidesign(x, y, column_design = NULL, cells = NULL, ...)
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

- cells:

  Optional logical matrix with dimensions identical to \`x\` and no
  missing values. \`NULL\` means no explicit mask; it is never inferred
  from missing values in \`x\`.

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

A multidesign object consists of three core components and an optional
mask: \* A data matrix where rows represent observations and columns
represent variables \* A design data frame containing experimental
factors and conditions for each observation \* Optional column metadata
describing properties of each variable \* An optional logical
cell-observation mask with the same dimensions as the data matrix

A multidesign object consists of three main components: \* A data matrix
where rows represent observations and columns represent variables \* A
design data frame containing experimental factors and conditions for
each observation \* Optional column metadata describing properties of
each variable \* An optional logical cell-observation mask aligned
exactly with the data matrix

The object maintains the relationship between these components while
providing methods for manipulation, subsetting, and analysis. A value of
\`NA\` in \`x\` is ordinary data unless the corresponding \`cells\`
entry is explicitly \`FALSE\`.

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

# Missing values do not imply missing cells; store a mask explicitly
cells <- matrix(TRUE, nrow(X), ncol(X))
cells[1, 1] <- FALSE
masked <- multidesign(X, Y, cells = cells)
has_cell_mask(masked)
#> [1] TRUE
cell_mask(masked)[1, 1]
#> [1] FALSE

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

# Store cell observation independently of data values
mask <- matrix(TRUE, nrow(X), ncol(X))
mask[1, 1] <- FALSE
masked_mds <- multidesign(X, Y, col_info, cells = mask)
cell_mask(masked_mds)[1:2, 1:2]
#>       [,1] [,2]
#> [1,] FALSE TRUE
#> [2,]  TRUE TRUE
```
