# Extract Design Information

Retrieves the design information (experimental metadata) associated with
observations in various object types. This function provides a
consistent interface for accessing design variables regardless of the
specific object structure.

## Usage

``` r
design(x, ...)

# S3 method for class 'multidesign'
design(x, ...)

# S3 method for class 'multiframe'
design(x, ...)
```

## Arguments

- x:

  The object containing design information

- ...:

  Additional arguments passed to methods: \* For hyperdesign objects:
  'block' parameter to specify which block's design to return

## Value

The design component of the object (typically a data frame or tibble),
or a list of designs for hyperdesign objects

## Details

The function behaves differently depending on the class of the input
object: \* For multidesign objects: Returns the design data frame
component \* For hyperdesign objects: Returns a list of design data
frames, one for each block (or a single design if a specific block is
requested) \* For multiframe objects: Returns the design tibble with
observation functions

Design information typically includes experimental factors, conditions,
subject IDs, or other metadata associated with each observation.

## See also

[`xdata`](https://bbuchsbaum.github.io/multidesign/reference/xdata.md)
for extracting the data matrix,
[`column_design`](https://bbuchsbaum.github.io/multidesign/reference/column_design.md)
for extracting column metadata

## Examples

``` r
# With a multidesign object
X <- matrix(rnorm(20*10), 20, 10)
Y <- data.frame(
  condition = rep(c("A", "B"), each=10),
  subject = rep(1:5, times=4)
)
mds <- multidesign(X, Y)
design_info <- design(mds)  # Returns the original data frame Y

# With a hyperdesign object
d1 <- multidesign(matrix(rnorm(10*5), 10, 5), 
                 data.frame(condition=rep(c("A","B"), 5), subject=1))
d2 <- multidesign(matrix(rnorm(10*5), 10, 5), 
                 data.frame(condition=rep(c("A","B"), 5), subject=2))
hd <- hyperdesign(list(d1, d2))
all_designs <- design(hd)  # Returns list of design data frames
block1_design <- design(hd, block=1)  # Returns just the first block's design
```
