# Create an Observation Object

Constructs a new lazy-evaluated observation object from various data
sources. An observation represents a single row, element, or vector from
a data source that can be accessed on demand.

## Usage

``` r
observation(x, i)

# S3 method for class 'deflist'
observation(x, i)

# Default S3 method
observation(x, i)

# S3 method for class 'matrix'
observation(x, i)

# S3 method for class 'list'
observation(x, i)
```

## Arguments

- x:

  The data source (matrix, list, vector, or other supported object)

- i:

  The index of the observation to extract

## Value

A function that, when called, returns the specified observation from the
data source

## See also

\[multiframe()\], \[obs_group()\]

## Examples

``` r
# From a matrix
X <- matrix(1:20, 5, 4)
obs1 <- observation(X, 2)  # Create observation for row 2
obs1()  # Retrieve the observation (returns row 2 as a 1xn matrix)
#>      [,1] [,2] [,3] [,4]
#> [1,]    2    7   12   17

# From a list
X_list <- list(a=1:5, b=6:10, c=11:15)
obs2 <- observation(X_list, 3)  # Create observation for 3rd element
obs2()  # Retrieve the observation (returns the 3rd element: 11:15)
#> [1] 11 12 13 14 15
```
