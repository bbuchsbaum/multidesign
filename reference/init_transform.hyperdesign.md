# Initialize Transformation for Hyperdesign

Method to initialize transformations (e.g., scaling, centering) for
hyperdesign objects. Each block in the hyperdesign gets its own
transformation object.

## Usage

``` r
# S3 method for class 'hyperdesign'
init_transform(x, X, ...)
```

## Arguments

- x:

  A hyperdesign object

- X:

  A preprocessing specification (e.g., from recipes package)

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
if (FALSE) { # \dontrun{
d1 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(cond = rep(c("A","B"), 5)))
d2 <- multidesign(matrix(rnorm(10*5), 10, 5),
                  data.frame(cond = rep(c("A","B"), 5)))
hd <- hyperdesign(list(d1, d2))
hd_transformed <- init_transform(hd, recipes::recipe(~ ., data = as.data.frame(d1$x)))
} # }
```
