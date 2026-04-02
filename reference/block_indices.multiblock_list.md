# Extract Block Indices from Multiblock Object

Retrieves the indices corresponding to a specific block in a multiblock
object. These indices can be used to access the corresponding rows or
columns in the combined matrix representation.

## Usage

``` r
# S3 method for class 'multiblock_list'
block_indices(x, i, ...)
```

## Arguments

- x:

  A multiblock object

- i:

  Integer specifying which block's indices to retrieve

- ...:

  Additional arguments (not used)

## Value

Integer vector of indices for the specified block

## See also

Other multiblock functions:
[`multiblock.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiblock.list.md),
[`t.multiblock_list()`](https://bbuchsbaum.github.io/multidesign/reference/t.multiblock_list.md)

## Examples

``` r
mb <- multiblock(list(matrix(1:10, 5, 2), matrix(11:20, 5, 2)))
block_indices(mb, 1)
#> [1] 1 2
block_indices(mb, 2)
#> [1] 3 4
```
