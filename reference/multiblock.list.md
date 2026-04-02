# Create a Multiblock Object from a List of Matrices

Creates a multiblock object from a list of matrices that share either
row or column dimensions. The function automatically determines whether
the matrices should be row-stacked or column-stacked based on their
shared dimensions.

## Usage

``` r
# S3 method for class 'list'
multiblock(x, ...)
```

## Arguments

- x:

  A list of matrices (base or Matrix class objects)

- ...:

  Additional arguments (not used)

## Value

A multiblock object with the following attributes: \* ind: matrix of
start/end indices for each block \* orient: orientation of stacking
("cstacked" or "rstacked") \* class: appropriate class labels for
dispatch

## Details

The function checks if the input matrices can be combined either by: \*
Row-stacking: all matrices must have the same number of columns \*
Column-stacking: all matrices must have the same number of rows

The resulting object maintains block structure information while
allowing operations across the entire combined matrix.

## See also

\[is_cstacked()\] for checking stacking orientation, \[block_indices()\]
for accessing block-specific indices

Other multiblock functions:
[`block_indices.multiblock_list()`](https://bbuchsbaum.github.io/multidesign/reference/block_indices.multiblock_list.md),
[`t.multiblock_list()`](https://bbuchsbaum.github.io/multidesign/reference/t.multiblock_list.md)

## Examples

``` r
# Create example matrices with shared row dimension (column-stacked)
X1 <- matrix(rnorm(20*3), 20, 3)
X2 <- matrix(rnorm(20*5), 20, 5)
X3 <- matrix(rnorm(20*4), 20, 4)
mb_c <- multiblock(list(X1, X2, X3))  # column-stacked
is_cstacked(mb_c)  # TRUE
#> [1] TRUE

# Create example matrices with shared column dimension (row-stacked)
Y1 <- matrix(rnorm(10*5), 10, 5)
Y2 <- matrix(rnorm(15*5), 15, 5)
Y3 <- matrix(rnorm(20*5), 20, 5)
mb_r <- multiblock(list(Y1, Y2, Y3))  # row-stacked
is_rstacked(mb_r)  # TRUE
#> [1] TRUE
```
