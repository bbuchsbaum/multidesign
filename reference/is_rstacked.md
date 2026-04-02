# Test if Multiblock Object is Row Stacked

Checks if a multiblock object is row-stacked, meaning that all component
matrices share the same number of columns. In row-stacked multiblock
objects, blocks are arranged one above another, with each block
potentially having a different number of rows.

## Usage

``` r
is_rstacked(x)

# S3 method for class 'multiblock_list'
is_rstacked(x)
```

## Arguments

- x:

  The multiblock object to test

## Value

Logical value: TRUE if the object is row-stacked, FALSE otherwise

## See also

[`is_cstacked`](https://bbuchsbaum.github.io/multidesign/reference/is_cstacked.md)
for testing if a multiblock object is column-stacked,
[`multiblock`](https://bbuchsbaum.github.io/multidesign/reference/multiblock.md)
for creating multiblock objects

## Examples

``` r
# Create row-stacked multiblock (matrices share column dimension)
Y1 <- matrix(rnorm(10*5), 10, 5)
Y2 <- matrix(rnorm(15*5), 15, 5)
mb <- multiblock(list(Y1, Y2))

# Test stacking orientation
is_rstacked(mb)  # Returns TRUE
#> [1] TRUE
is_cstacked(mb)  # Returns FALSE
#> [1] FALSE
```
