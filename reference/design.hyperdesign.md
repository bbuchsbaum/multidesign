# Extract Design Information from Hyperdesign

Retrieves design information from a hyperdesign object, either for all
blocks or for a specific block.

## Usage

``` r
# S3 method for class 'hyperdesign'
design(x, block, ...)
```

## Arguments

- x:

  A hyperdesign object

- block:

  Optional numeric index specifying which block's design to return

- ...:

  Additional arguments (not used)

## Value

If block is specified, returns design for that block; otherwise returns
a list of designs for all blocks

## See also

Other hyperdesign functions:
[`as_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/as_multidesign.md),
[`df_to_hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md),
[`hyperdesign.list()`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.list.md),
[`init_transform.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/init_transform.hyperdesign.md),
[`subset.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.hyperdesign.md),
[`xdata.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.hyperdesign.md)

## Examples

``` r
d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(y=1:10, condition=rep(c("A","B"), 5)))
d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(y=1:10, condition=rep(c("A","B"), 5)))
hd <- hyperdesign(list(d1, d2))

# Get all designs
all_designs <- design(hd)

# Get design for block 1
block1_design <- design(hd, block=1)
```
