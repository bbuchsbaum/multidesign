# Define a frame-native model design

A design specification preserves formulas and their environments. It
does not freeze a model matrix: coding, missing-value handling, and
multivariate block expansion occur when the specification is compiled
against a frame.

## Usage

``` r
design_spec(
  fixed,
  random = NULL,
  contrasts = NULL,
  na_action = c("fail", "omit")
)
```

## Arguments

- fixed:

  One-sided fixed-effects formula.

- random:

  Optional one-sided random-effects formula.

- contrasts:

  Optional named contrast list passed to \[stats::model.matrix()\].

- na_action:

  Missing-value policy, either \`"fail"\` or \`"omit"\`.

## Value

A serializable \`design_spec\`.
