# Create a Hyperdesign Object

Constructs a new \`hyperdesign\` object that encapsulates a collection
of \`multidesign\` objects. Used to model multi-group, multi-block, or
multi-view datasets, where each block/group/view is associated with a
matrix-variate response and an arbitrary design.

## Usage

``` r
hyperdesign(
  x,
  block_names = NULL,
  id = NULL,
  space = NULL,
  positional = FALSE,
  aggregate = NULL
)
```

## Arguments

- x:

  A list of multidesign objects. Each instance should represent a
  related block of data

- block_names:

  Optional character vector of names for each block

- id:

  Optional length-one character string naming the design column whose
  values identify corresponding entities across blocks.

- space:

  Optional column-space declaration: \`"common"\` when all blocks
  inhabit the same column space, or \`"block"\` when columns are
  block-specific.

- positional:

  Logical; if \`TRUE\`, declare positional row correspondence. This
  requires \`id = NULL\` and equal row counts in every block.

- aggregate:

  Optional duplicate-entity aggregation rule. \`NULL\` keeps duplicate
  IDs as an error; use \`"mean"\` or a scalar-returning function to
  aggregate within-block replicates explicitly.

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

## Details

A hyperdesign object represents a collection of related multivariate
datasets (multidesign instances) that share common design variables.
This structure is particularly useful for: \* Multiple subjects in an
experiment \* Multiple sessions or runs \* Multiple data modalities
(e.g., fMRI, EEG, behavioral) \* Multiple response measures

Rows retain the usual multidesign orientation: observations or
corresponding entities are rows, and measured variables are columns.
\`id\` gives a shared design column join semantics across blocks; it
does not transpose the stored matrices. \`common_vars\` describes shared
design-column names, whereas \`id\` identifies shared row values.
\`space = "common"\` additionally asserts that block columns represent
the same axes; \`space = "block"\` declares that they are
block-specific.

## See also

[`df_to_hyperdesign`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md)
for creating hyperdesign objects from data frames,
[`multidesign`](https://bbuchsbaum.github.io/multidesign/reference/multidesign.md)
for the underlying multidesign structure,
[`multiblock`](https://bbuchsbaum.github.io/multidesign/reference/multiblock.md)
for another multi-block data structure

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

# Duplicate entity IDs require an explicit aggregation rule
repeated <- multidesign(
  matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE),
  data.frame(entity = c("A", "A"))
)
aggregated <- hyperdesign(
  list(sample = repeated), id = "entity", space = "common",
  aggregate = "mean"
)
```
