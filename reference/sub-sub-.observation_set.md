# Extract Single Observation from Observation Set

Extracts a single observation from an observation_set using double
bracket indexing. The result is automatically evaluated.

## Usage

``` r
# S3 method for class 'observation_set'
x[[i]]
```

## Arguments

- x:

  An observation_set object

- i:

  Index of the observation to extract

## Value

The evaluated observation

## See also

Other multiframe functions:
[`[.observation_set()`](https://bbuchsbaum.github.io/multidesign/reference/sub-.observation_set.md),
[`multiframe.list()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.list.md),
[`multiframe.matrix()`](https://bbuchsbaum.github.io/multidesign/reference/multiframe.matrix.md),
[`obs_group()`](https://bbuchsbaum.github.io/multidesign/reference/obs_group.md),
[`split.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/split.multiframe.md),
[`subset.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/subset.multiframe.md),
[`summarize_by.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.multiframe.md),
[`xdata.multiframe()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.multiframe.md)
