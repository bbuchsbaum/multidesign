# Create a Multiframe from a List

Creates a multiframe object from a list of observations and associated
design information. A multiframe combines experimental design metadata
with lazy-evaluated observations. Each element in the list must have
consistent dimensions (same number of columns if matrices, or same
length if vectors).

## Usage

``` r
# S3 method for class 'list'
multiframe(x, y, ...)
```

## Arguments

- x:

  A list containing observations (matrices, arrays, or vectors)

- y:

  A data frame containing design variables. Must have same number of
  rows as length of x

- ...:

  Additional arguments (currently unused)

## Value

A multiframe object containing a design tibble with observation
functions

## See also

Other multiframe functions:
[`[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-.observation_set.md),
[`[[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-sub-.observation_set.md),
[`multiframe.matrix()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.matrix.md),
[`obs_group()`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md),
[`split.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/split.multiframe.md),
[`subset.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multiframe.md),
[`summarize_by.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multiframe.md),
[`xdata.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.multiframe.md)

## Examples

``` r
# Create list of observations (matrices with same number of columns)
x <- list(
  matrix(1:12, 3, 4),  # 3x4 matrix
  matrix(13:24, 3, 4), # 3x4 matrix
  matrix(25:36, 3, 4)  # 3x4 matrix
)

# Create design information
y <- data.frame(
  condition = c("A", "B", "C"),
  block = 1:3
)

# Create multiframe
mf <- multiframe(x, y)

# Access first observation
obs1 <- mf$design$.obs[[1]]()

# View design information
print(mf$design)
#> # A tibble: 3 × 4
#>   condition block .index .obs  
#>   <chr>     <int>  <int> <list>
#> 1 A             1      1 <fn>  
#> 2 B             2      2 <fn>  
#> 3 C             3      3 <fn>  
```
