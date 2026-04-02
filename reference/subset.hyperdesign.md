# Subset a Hyperdesign Object

Create a new hyperdesign object containing only the specified blocks.

## Usage

``` r
# S3 method for class 'hyperdesign'
subset(x, fexpr, ...)
```

## Arguments

- x:

  A hyperdesign object

- fexpr:

  Filter expression to apply to each block's design

- ...:

  Additional arguments (not used)

## Value

A new hyperdesign object containing only the selected blocks

## See also

Other hyperdesign functions:
[`as_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/as_multidesign.md),
[`design.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/design.hyperdesign.md),
[`df_to_hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md),
[`hyperdesign.list()`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.list.md),
[`init_transform.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/init_transform.hyperdesign.md),
[`xdata.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.hyperdesign.md)

## Examples

``` r
# Create example hyperdesign
d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(subject=1, condition=rep(c("A","B"), 5)))
d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(subject=2, condition=rep(c("A","B"), 5)))
d3 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(subject=3, condition=rep(c("A","B"), 5)))
hd <- hyperdesign(list(d1, d2, d3))

# Keep only condition A
subset_hd <- subset(hd, condition == "A")
```
