# Extract Column Design Information

Retrieves metadata associated with the columns (variables) of a
multidesign or related object. This function is particularly useful for
datasets where variables have associated metadata, such as feature
names, measurement types, or other variable-specific information.

## Usage

``` r
column_design(x, ...)

# S3 method for class 'hyperdesign'
column_design(x, block, ...)

# S3 method for class 'multidesign'
column_design(x, ...)
```

## Arguments

- x:

  The object containing column design information

- ...:

  Additional arguments passed to methods: \* For hyperdesign objects:
  'block' parameter to specify which block's column design to return

- block:

  Optional block index to get design for a specific block

## Value

The column design component of the object (typically a data frame or
tibble), or a list of column designs for hyperdesign objects

## Details

The function behaves differently depending on the class of the input
object: \* For multidesign objects: Returns the column_design data frame
component \* For hyperdesign objects: Returns a list of column design
data frames, one for each block (or a single column design if a specific
block is requested)

Column design information typically includes variable names, feature
types, regions of interest, or other metadata that describes the
variables rather than the observations.

## See also

[`design`](https://bbuchsbaum.github.io/multidesign/reference/design.md)
for extracting observation design information,
[`xdata`](https://bbuchsbaum.github.io/multidesign/reference/xdata.md)
for extracting the data matrix

## Examples

``` r
# With a multidesign object including column metadata
X <- matrix(rnorm(20*10), 20, 10)
Y <- data.frame(condition = rep(c("A", "B"), each=10))
col_info <- data.frame(
  feature = paste0("var", 1:10),
  type = rep(c("continuous", "categorical"), 5)
)
mds <- multidesign(X, Y, col_info)
col_metadata <- column_design(mds)

# With a hyperdesign object
d1 <- multidesign(matrix(rnorm(10*5), 10, 5), 
                 data.frame(condition=rep(c("A","B"), 5)),
                 data.frame(feature=paste0("var", 1:5)))
d2 <- multidesign(matrix(rnorm(10*5), 10, 5), 
                 data.frame(condition=rep(c("A","B"), 5)),
                 data.frame(feature=paste0("var", 1:5)))
hd <- hyperdesign(list(d1, d2))
all_col_designs <- column_design(hd)
block1_col_design <- column_design(hd, block=1)
```
