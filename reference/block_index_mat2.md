# Create Index Matrix for Block Structure

Internal function to compute start and end indices for blocks in a
multiblock object

## Usage

``` r
block_index_mat2(x, byrow = FALSE)
```

## Arguments

- x:

  List of matrices

- byrow:

  Logical; if TRUE compute row indices, if FALSE compute column indices

## Value

Matrix with two columns: 'start' and 'end' indices for each block
