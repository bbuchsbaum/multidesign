# Combine multiple multidesign objects

This function row-binds the observation matrices and design data frames
of several multidesign objects. All input multidesigns must share the
same column design. Optionally, an identifier column can be added to
track the source of each observation.

## Usage

``` r
bind_multidesign(..., .id = NULL)
```

## Arguments

- ...:

  multidesign objects to combine

- .id:

  optional name of an identifier column added to the design to indicate
  the origin multidesign

## Value

A new multidesign object containing the concatenated data

## Examples

``` r
X1 <- matrix(rnorm(10*5), 10, 5)
X2 <- matrix(rnorm(10*5), 10, 5)
md1 <- multidesign(X1, data.frame(cond = rep("A", 10)))
md2 <- multidesign(X2, data.frame(cond = rep("B", 10)))
combined <- bind_multidesign(md1, md2)
combined_id <- bind_multidesign(md1, md2, .id = "source")
```
