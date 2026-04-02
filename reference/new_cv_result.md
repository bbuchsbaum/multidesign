# Create a cv_result Object

Internal constructor for cv_result objects. Validates that scores is a
data frame with a \`.fold\` column and wraps it in a cv_result
structure.

## Usage

``` r
new_cv_result(scores, call = NULL)
```

## Arguments

- scores:

  A data frame containing cross-validation scores with a \`.fold\`
  column

- call:

  Optional call object recording how the result was created

## Value

A cv_result object
