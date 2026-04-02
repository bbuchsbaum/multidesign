# Construct a Hyperdesign Object

Creates a hyperdesign object, which represents a collection of related
multivariate datasets (multidesign instances) that share common design
variables. This class is particularly useful for modeling multi-block
data, where you want to analyze multiple related matrices, such as: \*
Multiple subjects in an experiment \* Multiple sessions or runs \*
Multiple data modalities (e.g., fMRI, EEG, behavioral) \* Multiple
response measures

## Usage

``` r
# S3 method for class 'list'
hyperdesign(x, block_names = NULL)
```

## Arguments

- x:

  A list of multidesign instances. Each instance should represent a
  related block of data

- block_names:

  Optional character vector of names for each block. If NULL, blocks
  will be automatically named as "block_1", "block_2", etc.

## Value

A hyperdesign object with the following components:

- blocks:

  List of multidesign objects

- block_names:

  Names of each block

- col_indices:

  Matrix of column start/end indices for each block

- row_indices:

  Matrix of row start/end indices for each block

## See also

[`df_to_hyperdesign`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md)
for creating hyperdesign objects from data frames,
[`multidesign`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md)
for the underlying multidesign structure,
[`design.hyperdesign`](https://bbuchsbaum.github.io/multidesign/reference/design.hyperdesign.md)
for extracting design information

Other hyperdesign functions:
[`as_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/as_multidesign.md),
[`design.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/design.hyperdesign.md),
[`df_to_hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md),
[`init_transform.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/init_transform.hyperdesign.md),
[`subset.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.hyperdesign.md),
[`xdata.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.hyperdesign.md)

## Examples

``` r
# Create three multidesign objects (e.g., for three subjects)
d1 <- multidesign(
  matrix(rnorm(10*20), 10, 20),
  data.frame(y=1:10, subject=1, run=rep(1:5, 2))
)
d2 <- multidesign(
  matrix(rnorm(10*20), 10, 20),
  data.frame(y=1:10, subject=2, run=rep(1:5, 2))
)
d3 <- multidesign(
  matrix(rnorm(10*20), 10, 20),
  data.frame(y=1:10, subject=3, run=rep(1:5, 2))
)

# Combine into a hyperdesign
hd <- hyperdesign(
  list(d1, d2, d3),
  block_names = c("subject1", "subject2", "subject3")
)
```
