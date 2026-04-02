# Create a Multiframe from a Matrix

Creates a multiframe object from a matrix of observations and associated
design information. Each row of the matrix becomes a lazy-evaluated
observation.

## Usage

``` r
# S3 method for class 'matrix'
multiframe(x, y, ...)
```

## Arguments

- x:

  A matrix where rows are observations and columns are variables

- y:

  A data frame containing design variables. Must have same number of
  rows as x

- ...:

  Additional arguments (currently unused)

## Value

A multiframe object containing a design tibble with observation
functions

## See also

Other multiframe functions:
[`[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-.observation_set.md),
[`[[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-sub-.observation_set.md),
[`multiframe.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.list.md),
[`obs_group()`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md),
[`split.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/split.multiframe.md),
[`subset.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multiframe.md),
[`summarize_by.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multiframe.md),
[`xdata.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.multiframe.md)

## Examples

``` r
# Create matrix of observations (10 observations x 4 variables)
x <- matrix(1:40, 10, 4)
colnames(x) <- paste0("var", 1:4)  # Optional: name the variables

# Create design information with two factors
y <- data.frame(
  condition = rep(c("A", "B"), each=5),
  subject = rep(1:5, times=2)
)

# Create multiframe
mf <- multiframe(x, y)  # Using the generic method

# Access first observation
first_obs <- mf$design$.obs[[1]]()
print(first_obs)  # Shows the first row of x as a matrix
#>      var1 var2 var3 var4
#> [1,]    1   11   21   31

# View design structure with observation functions
head(mf$design)
#> # A tibble: 6 × 4
#>   condition subject .index .obs  
#>   <chr>       <int>  <int> <list>
#> 1 A               1      1 <fn>  
#> 2 A               2      2 <fn>  
#> 3 A               3      3 <fn>  
#> 4 A               4      4 <fn>  
#> 5 A               5      5 <fn>  
#> 6 B               1      6 <fn>  
```
