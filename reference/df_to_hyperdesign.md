# Convert a Data Frame to a Hyperdesign Object

This function creates a hyperdesign object from a data frame by
splitting it according to a specified variable. It's particularly useful
when you have a wide-format data frame that needs to be converted into
multiple related multidesign objects, such as when dealing with multiple
subjects or sessions in an experiment.

## Usage

``` r
df_to_hyperdesign(data, design_vars, x_vars, split_var)
```

## Arguments

- data:

  A data frame or tibble containing both design variables and response
  variables

- design_vars:

  Character vector specifying the names of design variables (e.g.,
  conditions, factors)

- x_vars:

  Character vector specifying the names of response variables to extract

- split_var:

  Character string naming the variable to split the data on (e.g.,
  "subject" or "session")

## Value

A hyperdesign object containing multiple multidesign objects, one for
each unique value in split_var

## See also

[`hyperdesign`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.md)
for creating hyperdesign objects directly from multidesign objects,
[`multidesign`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md)
for the underlying multidesign structure,
[`design.hyperdesign`](https://bbuchsbaum.github.io/multidesign/reference/design.hyperdesign.md)
for extracting design information

Other hyperdesign functions:
[`as_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/as_multidesign.md),
[`design.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/design.hyperdesign.md),
[`hyperdesign.list()`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.list.md),
[`init_transform.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/init_transform.hyperdesign.md),
[`subset.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/subset.hyperdesign.md),
[`xdata.hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.hyperdesign.md)

## Examples

``` r
# Create a sample tibble with multiple subjects
sample_tibble <- tibble::tibble(
  felab = rep(1:2, each = 3),
  attention = rep(c("DA", "FA"), times = 3),
  basis = rep(c("basis01", "basis02", "basis03"), times = 2),
  subject = rep(1001:1002, each = 3),
  `1` = rnorm(6),  # response variable 1
  `2` = rnorm(6),  # response variable 2
  `3` = rnorm(6)   # response variable 3
)

# Convert to hyperdesign, splitting by subject
hd <- df_to_hyperdesign(
  data = sample_tibble,
  design_vars = c("felab", "attention", "basis"),
  x_vars = as.character(1:3),
  split_var = "subject"
)
```
