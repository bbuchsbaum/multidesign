# Test if Multiblock Object is Column Stacked

Checks if a multiblock object is column-stacked, meaning that all
component matrices share the same number of rows. In column-stacked
multiblock objects, blocks are arranged side by side, with each block
potentially having a different number of columns.

## Usage

``` r
is_cstacked(x)

# S3 method for class 'multiblock_list'
is_cstacked(x)
```

## Arguments

- x:

  The multiblock object to test

## Value

Logical value: TRUE if the object is column-stacked, FALSE otherwise

## See also

[`is_rstacked`](https://bbuchsbaum.github.io/multidesign/reference/is_rstacked.md)
for testing if a multiblock object is row-stacked,
[`multiblock`](https://bbuchsbaum.github.io/multidesign/reference/multiblock.md)
for creating multiblock objects

## Examples

``` r
# Create column-stacked multiblock (matrices share row dimension)
X1 <- matrix(rnorm(20*3), 20, 3)
X2 <- matrix(rnorm(20*5), 20, 5)
mb <- multiblock(list(X1, X2))

# Test stacking orientation
is_cstacked(mb)  # Returns TRUE
#> [1] TRUE
is_rstacked(mb)  # Returns FALSE
#> [1] FALSE
```
