# Create a Multiblock Object

Constructs a new multiblock object consisting of a set of stacked
submatrices sharing a row or column dimension. This structure is useful
for analyzing data with multiple related blocks of measurements while
preserving block structure information.

## Usage

``` r
multiblock(x, ...)
```

## Arguments

- x:

  A list of matrices sharing either row or column dimensions

- ...:

  Additional arguments passed to methods

## Value

A multiblock object with the following attributes:

- ind:

  Matrix of start/end indices for each block

- orient:

  Orientation of stacking ("cstacked" or "rstacked")

## Details

A multiblock object automatically determines whether matrices should be:
\* Row-stacked: all matrices must have the same number of columns \*
Column-stacked: all matrices must have the same number of rows

The resulting object maintains block structure information while
allowing operations across the entire combined matrix.

## See also

[`is_cstacked`](https://bbuchsbaum.github.io/multidesign/reference/is_cstacked.md)
for checking if a multiblock object is column-stacked,
[`is_rstacked`](https://bbuchsbaum.github.io/multidesign/reference/is_rstacked.md)
for checking if a multiblock object is row-stacked,
[`block_indices`](https://bbuchsbaum.github.io/multidesign/reference/block_indices.md)
for accessing block-specific indices

## Examples

``` r
# Create list of matrices with varying row dimensions (column-stacked)
X1 <- matrix(rnorm(20*3), 20, 3)
X2 <- matrix(rnorm(20*5), 20, 5)
X3 <- matrix(rnorm(20*4), 20, 4)
mb_c <- multiblock(list(X1, X2, X3))
is_cstacked(mb_c)  # TRUE
#> [1] TRUE

# Create list of matrices with varying column dimensions (row-stacked)
Y1 <- matrix(rnorm(10*5), 10, 5)
Y2 <- matrix(rnorm(15*5), 15, 5)
Y3 <- matrix(rnorm(20*5), 20, 5)
mb_r <- multiblock(list(Y1, Y2, Y3))
is_rstacked(mb_r)  # TRUE
#> [1] TRUE
```
