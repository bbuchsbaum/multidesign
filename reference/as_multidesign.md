# Convert to a Multidesign Object

Converts an object (e.g., a hyperdesign) into a single multidesign
object by collapsing its structure.

Converts a hyperdesign object into a single multidesign by row-stacking
the data matrices and combining the design data frames. All blocks must
have the same number of columns and identical column designs.

## Usage

``` r
as_multidesign(x, ...)

# S3 method for class 'hyperdesign'
as_multidesign(x, .id = NULL, ...)
```

## Arguments

- x:

  A hyperdesign object

- ...:

  Additional arguments (not used)

- .id:

  Optional character string. If provided, a column with this name is
  added to the design to identify the source block.

## Value

A multidesign object

A single multidesign object

## See also

[`hyperdesign`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.md)
for creating hyperdesign objects,
[`multidesign`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md)
for the multidesign class

Other hyperdesign functions:
[`design.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/design.hyperdesign.md),
[`df_to_hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md),
[`hyperdesign.list()`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.list.md),
[`init_transform.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/init_transform.hyperdesign.md),
[`subset.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.hyperdesign.md),
[`xdata.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.hyperdesign.md)

## Examples

``` r
d1 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(condition = rep(c("A","B"), 5)))
d2 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(condition = rep(c("A","B"), 5)))
hd <- hyperdesign(list(d1, d2))

# Collapse hyperdesign to a single multidesign
md <- as_multidesign(hd)

d1 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(condition = rep(c("A","B"), 5)))
d2 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(condition = rep(c("A","B"), 5)))
hd <- hyperdesign(list(d1, d2))
md <- as_multidesign(hd)
```
