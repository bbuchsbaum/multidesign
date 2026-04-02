# Extract Multiple Observations from Observation Set

Extracts multiple observations from an observation_set using single
bracket indexing. All selected observations are automatically evaluated.

## Usage

``` r
# S3 method for class 'observation_set'
x[i]
```

## Arguments

- x:

  An observation_set object

- i:

  Indices of observations to extract

## Value

List of evaluated observations

## See also

Other multiframe functions:
[`[[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-sub-.observation_set.md),
[`multiframe.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.list.md),
[`multiframe.matrix()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.matrix.md),
[`obs_group()`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md),
[`split.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/split.multiframe.md),
[`subset.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multiframe.md),
[`summarize_by.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multiframe.md),
[`xdata.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.multiframe.md)
