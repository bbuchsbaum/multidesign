# Create a Group of Observations from Matrix or List Data

Creates an observation group object that provides a unified interface
for accessing individual observations from different data structures
(matrices, lists, or deflists). Each observation is wrapped in a
function that provides lazy evaluation.

## Usage

``` r
obs_group(X, fun = NULL, ind = NULL)
```

## Arguments

- X:

  A matrix, list, or deflist containing the data

- fun:

  Optional function to apply to each observation

- ind:

  Optional vector of indices for the observations. If NULL, uses
  sequential indices

## Value

An observation_set object containing lazy-evaluated observations

## Details

Create an Observation Group

## See also

Other multiframe functions:
[`[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-.observation_set.md),
[`[[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-sub-.observation_set.md),
[`multiframe.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.list.md),
[`multiframe.matrix()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.matrix.md),
[`split.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/split.multiframe.md),
[`subset.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multiframe.md),
[`summarize_by.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multiframe.md),
[`xdata.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.multiframe.md)

## Examples

``` r
# From a matrix
X <- matrix(rnorm(100), 20, 5)
obs <- obs_group(X)

# From a list
X_list <- list(a=1:5, b=6:10, c=11:15)
obs_list <- obs_group(X_list)
```
