# Transpose a Multiblock Object

Transposes each matrix in a multiblock object and returns a new
multiblock object. The stacking orientation will be switched
(row-stacked becomes column-stacked and vice versa).

## Usage

``` r
# S3 method for class 'multiblock_list'
t(x)
```

## Arguments

- x:

  A multiblock object

## Value

A new multiblock object with transposed matrices

## See also

Other multiblock functions:
[`block_indices.multiblock_list()`](https://bbuchsbaum.github.io/multidesign/reference/block_indices.multiblock_list.md),
[`multiblock.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiblock.list.md)

## Examples

``` r
mb <- multiblock(list(matrix(1:10, 5, 2), matrix(11:20, 5, 2)))
mb_t <- t(mb)
is_cstacked(mb_t)
#> [1] TRUE
is_rstacked(mb_t)
#> [1] FALSE
```
