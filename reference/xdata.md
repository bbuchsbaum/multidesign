# Extract Data Matrix

Retrieves the data matrix component from various object types in the
package. This function provides a consistent interface for accessing the
underlying data regardless of the specific object structure.

## Usage

``` r
xdata(x, ...)

# S3 method for class 'multidesign'
xdata(x, ...)
```

## Arguments

- x:

  The object containing the data matrix

- ...:

  Additional arguments passed to methods: \* For hyperdesign objects:
  'block' parameter to specify which block's data to return

## Value

The data matrix component of the object, or a list of matrices for
hyperdesign objects

## Details

The function behaves differently depending on the class of the input
object: \* For multidesign objects: Returns the data matrix component \*
For hyperdesign objects: Returns a list of data matrices, one for each
block (or a single matrix if a specific block is requested) \* For
multiframe objects: Returns the combined data from all observations

## See also

[`design`](https://bbuchsbaum.github.io/multidesign/reference/design.md)
for extracting design information,
[`column_design`](https://bbuchsbaum.github.io/multidesign/reference/column_design.md)
for extracting column metadata

## Examples

``` r
# With a multidesign object
X <- matrix(rnorm(20*10), 20, 10)
Y <- data.frame(group = rep(letters[1:4], each=5))
mds <- multidesign(X, Y)
X_data <- xdata(mds)  # Returns the original matrix X

# With a hyperdesign object
d1 <- multidesign(matrix(rnorm(10*5), 10, 5), data.frame(subject=rep(1,10)))
d2 <- multidesign(matrix(rnorm(10*5), 10, 5), data.frame(subject=rep(2,10)))
hd <- hyperdesign(list(d1, d2))
all_data <- xdata(hd)  # Returns list of matrices
block1_data <- xdata(hd, block=1)  # Returns just the first block's matrix
```
