# Initialize Transformation for Hyperdesign

Method to initialize transformations (e.g., scaling, centering) for
hyperdesign objects. Each block in the hyperdesign gets its own
transformation object. Hyperdesigns with partial cell masks are rejected
because preprocessors do not declare how masks transform. All-\`TRUE\`
masks are reshaped with the output, while absent masks remain absent.

## Usage

``` r
# S3 method for class 'hyperdesign'
init_transform(x, X, ...)
```

## Arguments

- x:

  A hyperdesign object

- X:

  A preprocessing specification supported by \`multivarious\`, such as
  \`multivarious::center()\`.

- ...:

  Additional arguments (not used)

## Value

A hyperdesign with transformed data and a `preproc` attribute containing
the fitted preprocessing objects

## See also

Other hyperdesign functions:
[`as_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/as_multidesign.md),
[`design.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/design.hyperdesign.md),
[`df_to_hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md),
[`hyperdesign.list()`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.list.md),
[`subset.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.hyperdesign.md),
[`xdata.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.hyperdesign.md)

## Examples

``` r
d1 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(cond = rep(c("A","B"), 5)))
d2 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(cond = rep(c("A","B"), 5)))
hd <- hyperdesign(list(d1, d2))
hd_transformed <- init_transform(hd, multivarious::center())
```
