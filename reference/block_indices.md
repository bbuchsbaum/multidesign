# Get Block Indices from a Multiblock Object

This is a re-export of the block_indices function from the multivarious
package. See
[`multivarious::block_indices`](https://bbuchsbaum.github.io/multivarious/reference/block_indices.html)
for full documentation.

## Usage

``` r
block_indices(x, ...)

# S3 method for class 'hyperdesign'
block_indices(x, i, byrow = FALSE, ...)
```

## Arguments

- x:

  The object to get block indices from

- ...:

  Additional arguments passed to methods

- i:

  the block number

- byrow:

  if true, return row-oriented indices

## Value

An integer vector or list of integer vectors with indices for the
requested block(s)

## Details

Get Block Indices

## Examples

``` r
mb <- multiblock(list(matrix(1:10, 5, 2), matrix(11:20, 5, 2)))
block_indices(mb, 1)
#> [1] 1 2
block_indices(mb, 2)
#> [1] 3 4
```
