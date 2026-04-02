# Extract Data Matrix from Hyperdesign

Get the data matrix component from a hyperdesign object, either for all
blocks or for a specific block.

## Usage

``` r
# S3 method for class 'hyperdesign'
xdata(x, block, ...)
```

## Arguments

- x:

  A hyperdesign object

- block:

  Optional numeric index specifying which block's data to return

- ...:

  Additional arguments (not used)

## Value

If block is specified, returns data matrix for that block; otherwise
returns a list of data matrices for all blocks

## See also

\[design.hyperdesign()\], \[column_design.hyperdesign()\]

Other hyperdesign functions:
[`as_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/as_multidesign.md),
[`design.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/design.hyperdesign.md),
[`df_to_hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md),
[`hyperdesign.list()`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.list.md),
[`init_transform.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/init_transform.hyperdesign.md),
[`subset.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.hyperdesign.md)

## Examples

``` r
# Create example hyperdesign
d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(y=1:10, condition=rep(c("A","B"), 5)))
d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(y=1:10, condition=rep(c("A","B"), 5)))
hd <- hyperdesign(list(d1, d2))

# Get data from all blocks
all_data <- xdata(hd)

# Get data from block 1
block1_data <- xdata(hd, block=1)
```
