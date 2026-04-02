# Class Reference and Design Patterns

## Overview

The `multidesign` package provides a family of S3 classes for
representing multivariate experimental data. This vignette serves as a
reference guide, documenting each class, its purpose, key methods, and
typical usage patterns.

### Class hierarchy

    multidesign ─────────────────────── Core class: matrix + row design + column design
         │
         ├── reduced_multidesign ────── After dimensionality reduction
         │
         └── hyperdesign ────────────── Collection of multidesign blocks
                  │
                  └── as_multidesign() → Collapse back to single multidesign

    multiframe ──────────────────────── Lazy-evaluated observations + design

    multiblock ──────────────────────── Stacked matrices (row or column oriented)

    foldlist ────────────────────────── Cross-validation fold container

    cv_result ───────────────────────── Cross-validation results

### Choosing the right class

    Do you have multiple subjects/sessions/modalities?
    ├── YES → hyperdesign
    └── NO → Is your data expensive to load?
             ├── YES → multiframe
             └── NO → multidesign

    Do you need to stack matrices without design info?
    └── YES → multiblock

## multidesign

### Purpose

A `multidesign` is the fundamental class in the package. It binds
together:

- **Data matrix** (`x`): rows are observations, columns are variables
- **Row design** (`design`): a tibble describing each observation
- **Column design** (`column_design`): a tibble describing each variable
  (optional)

### When to use

- Standard single-block multivariate analysis
- You need to subset or split data by experimental factors
- You want cross-validation that respects design structure
- Your variables have metadata (ROI names, feature types, etc.)

### When NOT to use

- You have multiple subjects/sessions → use `hyperdesign`
- Your data is too large to load at once → use `multiframe`
- You just need to stack matrices without design info → use `multiblock`

### Construction

``` r
# Data matrix: 30 observations x 10 variables
X <- matrix(rnorm(30 * 10), 30, 10)

# Row design: what describes each observation
row_design <- tibble(
  condition = rep(c("A", "B", "C"), each = 10),
  block = rep(1:3, 10)
)

# Column design: what describes each variable (optional)
col_design <- tibble(
  variable = paste0("V", 1:10),
  group = rep(c("group1", "group2"), 5)
)

# Create multidesign
md <- multidesign(X, row_design, col_design)
md
#> 
#> === Multidesign Object ===
#> 
#> Data Matrix: 
#>    30 observations x 10 variables 
#> 
#> Design Variables: 
#>   * condition: 3 levels (A, B, C)
#>   * block: 3 levels (1, 2, 3)
#> 
#> Column Metadata:
#>   * variable: 10 levels (V1, V2, V3...V9, V10)
#>   * group: 2 levels (group1, group2)
#> 
#> =======================
#> 
```

### Key methods

| Method                     | Description                                    |
|----------------------------|------------------------------------------------|
| `xdata(x)`                 | Extract the data matrix                        |
| `design(x)`                | Extract row design (without internal `.index`) |
| `column_design(x)`         | Extract column design                          |
| `subset(x, expr)`          | Filter rows by design expression               |
| `split(x, ...)`            | Split into list by design variables            |
| `select_variables(x, ...)` | Filter columns by column_design expression     |
| `summarize_by(x, ...)`     | Compute group summaries                        |
| `fold_over(x, ...)`        | Create cross-validation folds                  |
| `cv_rows(x, rows)`         | Create folds from explicit row indices         |
| `split_indices(x, ...)`    | Get row indices grouped by design variables    |

### Examples

``` r
# Extract components
dim(xdata(md))
#> [1] 30 10
head(design(md), 3)
#> # A tibble: 3 × 2
#>   condition block
#>   <chr>     <int>
#> 1 A             1
#> 2 A             2
#> 3 A             3
column_design(md)
#> # A tibble: 10 × 2
#>    variable group 
#>    <chr>    <chr> 
#>  1 V1       group1
#>  2 V2       group2
#>  3 V3       group1
#>  4 V4       group2
#>  5 V5       group1
#>  6 V6       group2
#>  7 V7       group1
#>  8 V8       group2
#>  9 V9       group1
#> 10 V10      group2

# Subset by design
md_A <- subset(md, condition == "A")
nrow(xdata(md_A))
#> [1] 10

# Select variables by column design
md_g1 <- select_variables(md, group == "group1")
ncol(xdata(md_g1))
#> [1] 5

# Summarize by condition
md_means <- summarize_by(md, condition)
design(md_means)
#> # A tibble: 3 × 1
#>   condition
#>   <chr>    
#> 1 A        
#> 2 B        
#> 3 C
dim(xdata(md_means))
#> [1]  3 10

# Create CV folds by condition
folds <- fold_over(md, condition)
length(folds)
#> [1] 3

# Preserve original source row ids when needed downstream
folds_with_ids <- fold_over(md, condition, preserve_row_ids = TRUE)
head(folds_with_ids[[1]]$assessment$design$.orig_index)
#> [1] 1 2 3 4 5 6
```

------------------------------------------------------------------------

## reduced_multidesign

### Purpose

A `reduced_multidesign` extends `multidesign` to store the result of
dimensionality reduction (e.g., PCA). It contains the reduced data and
the projector used for the transformation.

### When to use

- After applying `reduce()` to a multidesign
- When you need to project new data using the same transformation

### Construction

``` r
# Note: requires multivarious package for full functionality
# Basic structure shown here
md_small <- multidesign(
  matrix(rnorm(50 * 20), 50, 20),
  tibble(group = rep(letters[1:5], 10))
)

# reduce() performs PCA by default
# rd <- reduce(md_small, nc = 3)
# rd
```

The `reduce()` function requires the `multivarious` package. When
called, it returns a `reduced_multidesign` with:

- Reduced data matrix (observations × components)
- Original design information
- The projector object for transforming new data

------------------------------------------------------------------------

## hyperdesign

### Purpose

A `hyperdesign` manages multiple `multidesign` objects (called “blocks”)
that share common design variables. This is the natural structure for
multi-subject, multi-session, or multi-modal experiments.

### When to use

- Multiple subjects, each with their own data matrix
- Multiple sessions or runs that should be analyzed together
- Multi-view data (e.g., fMRI + EEG from the same experiment)
- You want leave-one-block-out cross-validation

### When NOT to use

- Single-subject analysis → use `multidesign`
- Blocks don’t share any design structure → consider separate
  `multidesign` objects

### Construction

``` r
# Create individual subject data
subj1 <- multidesign(
  matrix(rnorm(20 * 5), 20, 5),
  tibble(condition = rep(c("X", "Y"), 10), trial = 1:20)
)

subj2 <- multidesign(
  matrix(rnorm(20 * 5), 20, 5),
  tibble(condition = rep(c("X", "Y"), 10), trial = 1:20)
)

subj3 <- multidesign(
  matrix(rnorm(20 * 5), 20, 5),
  tibble(condition = rep(c("X", "Y"), 10), trial = 1:20)
)

# Combine into hyperdesign
hd <- hyperdesign(
  list(subj1, subj2, subj3),
  block_names = c("subject_1", "subject_2", "subject_3")
)
hd
#> 
#> === Hyperdesign Object ===
#> 
#> Number of blocks:  3 
#> 
#> +- Block  1  (subject_1)  -----------------
#> |  Dimensions: 20 x 5 
#> |  Design Variables: condition, trial 
#> |  Design Structure: 
#> |   * condition: 2 levels (X, Y)
#> |   * trial: 20 levels (1, 2, 3...19, 20)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> +- Block  2  (subject_2)  -----------------
#> |  Dimensions: 20 x 5 
#> |  Design Variables: condition, trial 
#> |  Design Structure: 
#> |   * condition: 2 levels (X, Y)
#> |   * trial: 20 levels (1, 2, 3...19, 20)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> +- Block  3  (subject_3)  -----------------
#> |  Dimensions: 20 x 5 
#> |  Design Variables: condition, trial 
#> |  Design Structure: 
#> |   * condition: 2 levels (X, Y)
#> |   * trial: 20 levels (1, 2, 3...19, 20)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> =======================
#> 
```

### Key methods

| Method                             | Description                                   |
|------------------------------------|-----------------------------------------------|
| `xdata(x)` / `xdata(x, block=i)`   | Extract data (all blocks or specific block)   |
| `design(x)` / `design(x, block=i)` | Extract design (all blocks or specific block) |
| `column_design(x)`                 | Extract column design                         |
| `subset(x, expr)`                  | Filter all blocks by expression               |
| `select_variables(x, ...)`         | Filter columns in all blocks                  |
| `fold_over(x)`                     | Leave-one-block-out folds                     |
| `fold_over(x, var)`                | Within-block folds by variable                |
| `cv_rows(x, rows)`                 | Create synchronized row-index folds           |
| `as_multidesign(x)`                | Collapse to single multidesign                |
| `block_indices(x, i)`              | Get indices for block i                       |

### Examples

``` r
# Access individual blocks
length(hd)
#> [1] 3
xdata(hd, block = 1)[1:3, ]
#>            [,1]       [,2]       [,3]      [,4]      [,5]
#> [1,] -0.7497258 -0.5836394 -0.2185938 1.9001363 0.4463130
#> [2,] -0.3216061 -1.9940788  1.4599659 0.1100091 0.4218847
#> [3,] -1.1477708  1.9022098 -0.5820599 1.1403868 0.4424648

# Access all designs
designs <- design(hd)
length(designs)
#> [1] 3

# Leave-one-block-out cross-validation
loso <- fold_over(hd)
length(loso)
#> [1] 3

# Each fold: analysis is a hyperdesign (n-1 blocks), assessment is a multidesign (1 block)
f1 <- loso[[1]]
length(f1$analysis)
#> [1] 2
class(f1$assessment)
#> [1] "multidesign"

# Within-block CV by condition
within_folds <- fold_over(hd, condition)
length(within_folds)
#> [1] 6

# Explicit synchronized row holdouts
row_folds <- cv_rows(
  hd,
  rows = list(list(subject_1 = 1:2, subject_2 = 1:2)),
  preserve_row_ids = TRUE
)
row_folds[[1]]$assessment[[1]]$design$.orig_index
#> [1] 1 2

# Collapse to single multidesign
md_collapsed <- as_multidesign(hd, .id = "subject")
nrow(xdata(md_collapsed))
#> [1] 60
table(design(md_collapsed)$subject)
#> 
#> subject_1 subject_2 subject_3 
#>        20        20        20
```

### Alternative construction: from data frame

``` r
# If your data is in a wide-format data frame
wide_df <- tibble(
  subject = rep(1:2, each = 5),
  condition = rep(c("A", "B", "A", "B", "A"), 2),
  v1 = rnorm(10),
  v2 = rnorm(10),
  v3 = rnorm(10)
)

hd_from_df <- df_to_hyperdesign(
  wide_df,
  design_vars = c("condition"),
  x_vars = c("v1", "v2", "v3"),
  split_var = "subject"
)
hd_from_df
#> 
#> === Hyperdesign Object ===
#> 
#> Number of blocks:  2 
#> 
#> +- Block  1  (1)  -----------------
#> |  Dimensions: 5 x 3 
#> |  Design Variables: condition 
#> |  Design Structure: 
#> |   * condition: 2 levels (A, B)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> +- Block  2  (2)  -----------------
#> |  Dimensions: 5 x 3 
#> |  Design Variables: condition 
#> |  Design Structure: 
#> |   * condition: 2 levels (A, B)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> =======================
#> 
```

------------------------------------------------------------------------

## multiframe

### Purpose

A `multiframe` is similar to a `multidesign` but uses lazy evaluation
for observations. Each observation is stored as a function that returns
the data when called. This is useful when data is expensive to load or
too large to hold entirely in memory.

### When to use

- Data files are large and you only need subsets at a time
- Observations are loaded from disk on demand
- You want to manipulate design metadata without loading all data

### When NOT to use

- Data fits easily in memory → use `multidesign` (simpler)
- You need to operate on the full matrix frequently → use `multidesign`

### Construction

``` r
# From a matrix (each row becomes a lazy observation)
X <- matrix(1:40, 10, 4)
design_info <- tibble(
  group = rep(c("A", "B"), 5),
  id = 1:10
)

mf <- multiframe(X, design_info)
mf
#> 
#> === Multiframe Object ===
#> 
#>  Number of Observations: 10 
#> 
#> Design Variables: 
#>   - group: 2 levels (A, B)
#>   - id: 10 levels (1, 2, 3...9, 10)
#> 
#> Sample Observations: 
#>   - Observation 2: 1 x 4 matrix
#>   - Observation 6: 1 x 4 matrix
#>   - Observation 10: 1 x 4 matrix
#>    -   ... and 7 more observations 
#> 
#> ===================
#> 

# From a list (each element is a lazy observation)
obs_list <- list(
  matrix(1:12, 3, 4),
  matrix(13:24, 3, 4),
  matrix(25:36, 3, 4)
)
mf_list <- multiframe(obs_list, tibble(condition = c("X", "Y", "Z")))
```

### Key methods

| Method                  | Description                                    |
|-------------------------|------------------------------------------------|
| `design(x)`             | Extract design (without `.obs` and `.index`)   |
| `xdata(x)`              | Materialize all observations into a matrix     |
| `subset(x, expr)`       | Filter by design (lazy observations preserved) |
| `split(x, ...)`         | Split by design variables                      |
| `split_indices(x, ...)` | Get indices grouped by design                  |
| `summarize_by(x, ...)`  | Compute group summaries                        |
| `fold_over(x, ...)`     | Create cross-validation folds                  |
| `cv_rows(x, rows)`      | Create folds from explicit row indices         |

### Examples

``` r
# Design without internal columns
design(mf)
#> # A tibble: 10 × 2
#>    group    id
#>    <chr> <int>
#>  1 A         1
#>  2 B         2
#>  3 A         3
#>  4 B         4
#>  5 A         5
#>  6 B         6
#>  7 A         7
#>  8 B         8
#>  9 A         9
#> 10 B        10

# Materialize all data (evaluates all lazy observations)
mat <- xdata(mf)
dim(mat)
#> [1] 10  4

# Subset (lazy observations remain lazy)
mf_A <- subset(mf, group == "A")
nrow(mf_A$design)
#> [1] 5

# Access individual observations lazily
mf$design$.obs[[1]]()
#>      [,1] [,2] [,3] [,4]
#> [1,]    1   11   21   31

# Create CV folds
mf_folds <- fold_over(mf, group)
length(mf_folds)
#> [1] 2

# Preserve original source row ids in fold outputs
mf_folds_with_ids <- fold_over(mf, group, preserve_row_ids = TRUE)
mf_folds_with_ids[[1]]$assessment$design$.orig_index
#> [1] 1 3 5 7 9
```

------------------------------------------------------------------------

## multiblock

### Purpose

A `multiblock` represents a collection of matrices that share either a
row dimension (column-stacked) or a column dimension (row-stacked). It’s
a lower-level structure than `hyperdesign`, useful when you need to work
with the combined matrix directly.

### When to use

- Matrices that should be stacked into one large matrix
- Block-wise operations on combined data
- Input to multi-block statistical methods

### When NOT to use

- You need design metadata → use `hyperdesign`
- Matrices don’t share a dimension → can’t be stacked

### Construction

``` r
# Column-stacked: same number of rows, different columns
A1 <- matrix(1:12, 4, 3)
A2 <- matrix(13:20, 4, 2)
mb_col <- multiblock(list(A1, A2))
mb_col
#> 
#> === Multiblock Object ===
#> 
#> Number of blocks:  2 
#> Orientation:  cstacked 
#> 
#> Shared dimension:  4 rows 
#> Block-specific columns:  3, 2 
#> Total columns:  5 
#> 
#> =======================
#> 

# Row-stacked: same number of columns, different rows
B1 <- matrix(1:10, 2, 5)
B2 <- matrix(11:25, 3, 5)
mb_row <- multiblock(list(B1, B2))
mb_row
#> 
#> === Multiblock Object ===
#> 
#> Number of blocks:  2 
#> Orientation:  rstacked 
#> 
#> Shared dimension:  5 columns 
#> Block-specific rows:  2, 3 
#> Total rows:  5 
#> 
#> =======================
#> 
```

### Key methods

| Method                | Description                      |
|-----------------------|----------------------------------|
| `is_cstacked(x)`      | Check if column-stacked          |
| `is_rstacked(x)`      | Check if row-stacked             |
| `block_indices(x, i)` | Get indices for block i          |
| `t(x)`                | Transpose (switches orientation) |

### Examples

``` r
# Check orientation
is_cstacked(mb_col)
#> [1] TRUE
is_rstacked(mb_row)
#> [1] TRUE

# Get block indices
block_indices(mb_col, 1)
#> [1] 1 2 3
block_indices(mb_col, 2)
#> [1] 4 5

# Transpose switches orientation
mb_t <- t(mb_col)
is_rstacked(mb_t)
#> [1] TRUE
```

------------------------------------------------------------------------

## foldlist

### Purpose

A `foldlist` is returned by
[`fold_over()`](https://bbuchsbaum.github.io/multidesign/reference/fold_over.md)
or
[`cv_rows()`](https://bbuchsbaum.github.io/multidesign/reference/cv_rows.md)
and contains cross-validation folds. Each fold has an analysis set
(training data) and an assessment set (test data). \## Structure

Each element of a foldlist contains:

- `analysis`: The training data (same class as input)
- `assessment`: The test data (same class as input, or `multidesign` for
  hyperdesign assessment)
- `held_out`: Information about what was held out (optional)

When folds are created with `preserve_row_ids = TRUE`, the fold outputs
also carry a reserved `.orig_index` design column, and
`held_out$row_ids` records the same source-row identities explicitly.

### Examples

``` r
md <- multidesign(
  matrix(rnorm(40), 10, 4),
  tibble(group = rep(c("A", "B"), 5))
)

folds <- fold_over(md, group)
folds
#> 
#> === Cross-Validation Folds ===
#> 
#>  Number of Folds: 2 
#> 
#> Fold 1: 
#>   * Analysis Set: 5 observations
#>   * Assessment Set: 5 observations
#> 
#> Fold 2: 
#>   * Analysis Set: 5 observations
#>   * Assessment Set: 5 observations
#> 
#> =============================
#> 

# Access a specific fold
f1 <- folds[[1]]
names(f1)
#> [1] "analysis"   "assessment"

# Training and test sets
nrow(xdata(f1$analysis))
#> [1] 5
nrow(xdata(f1$assessment))
#> [1] 5

# The foldframe attribute contains fold metadata
attr(folds, "foldframe")
#> # A tibble: 2 × 4
#>   group indices   .splitvar .fold
#>   <chr> <list>    <chr>     <int>
#> 1 A     <int [5]> A             1
#> 2 B     <int [5]> B             2

# Explicit row folds work too
row_folds <- cv_rows(md, rows = list(1:2), preserve_row_ids = TRUE)
row_folds[[1]]$held_out$row_ids
#> [1] 1 2
```

------------------------------------------------------------------------

## cv_result

### Purpose

A `cv_result` object stores the results of cross-validation executed via
[`cross_validate()`](https://bbuchsbaum.github.io/multidesign/reference/cross_validate.md).
It contains per-fold scores and provides methods for summarizing
performance.

### Structure

- `scores`: A tibble with a `.fold` column and one or more metric
  columns
- `call`: The original call to
  [`cross_validate()`](https://bbuchsbaum.github.io/multidesign/reference/cross_validate.md)
  (optional)

### Key methods

| Method       | Description                                          |
|--------------|------------------------------------------------------|
| `print(x)`   | Display fold count and metric summaries              |
| `summary(x)` | Return tibble with mean/sd/min/max/median per metric |

### Examples

``` r
# Run cross-validation
md <- multidesign(
  matrix(rnorm(60), 15, 4),
  tibble(condition = rep(c("A", "B", "C"), 5))
)
folds <- fold_over(md, condition)

cv_res <- cross_validate(
  folds,
  fit_fn = function(analysis) {
    list(center = colMeans(xdata(analysis)))
  },
  score_fn = function(model, assessment) {
    mat <- xdata(assessment)
    pred <- matrix(model$center, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
    c(mse = mean((mat - pred)^2), mae = mean(abs(mat - pred)))
  }
)

# Print summary
cv_res
#> 
#> === Cross-Validation Result ===
#> 
#> Folds: 3 
#> 
#> Metrics:
#>    mse : mean = 0.6518 , sd = 0.1936 
#>    mae : mean = 0.6409 , sd = 0.0508 
#> 
#> ===============================
#> 

# Detailed summary
summary(cv_res)
#> # A tibble: 2 × 6
#>   metric  mean     sd   min   max median
#>   <chr>  <dbl>  <dbl> <dbl> <dbl>  <dbl>
#> 1 mse    0.652 0.194  0.520 0.874  0.562
#> 2 mae    0.641 0.0508 0.598 0.697  0.627

# Access raw scores
cv_res$scores
#> # A tibble: 3 × 3
#>     mse   mae .fold
#>   <dbl> <dbl> <int>
#> 1 0.562 0.598     1
#> 2 0.520 0.627     2
#> 3 0.874 0.697     3
```

------------------------------------------------------------------------

## Combining and converting between classes

### bind_multidesign: Combine multiple multidesign objects

``` r
md1 <- multidesign(matrix(1:10, 5, 2), tibble(g = rep("A", 5)))
md2 <- multidesign(matrix(11:20, 5, 2), tibble(g = rep("B", 5)))

# Bind with source tracking
combined <- bind_multidesign(md1, md2, .id = "source")
design(combined)
#> # A tibble: 10 × 2
#>    g     source
#>    <chr>  <int>
#>  1 A          1
#>  2 A          1
#>  3 A          1
#>  4 A          1
#>  5 A          1
#>  6 B          2
#>  7 B          2
#>  8 B          2
#>  9 B          2
#> 10 B          2

# Also accepts hyperdesign (expands blocks first)
hd <- hyperdesign(list(md1, md2))
md3 <- multidesign(matrix(21:30, 5, 2), tibble(g = rep("C", 5)))
combined2 <- bind_multidesign(hd, md3)
nrow(xdata(combined2))
#> [1] 15
```

### as_multidesign: Collapse hyperdesign

``` r
hd <- hyperdesign(
  list(
    multidesign(matrix(1:10, 5, 2), tibble(cond = c("X", "X", "Y", "Y", "Y"))),
    multidesign(matrix(11:20, 5, 2), tibble(cond = c("X", "X", "Y", "Y", "Y")))
  ),
  block_names = c("block1", "block2")
)

# Collapse to single multidesign
md_flat <- as_multidesign(hd, .id = "block")
design(md_flat)
#> # A tibble: 10 × 2
#>    cond  block 
#>    <chr> <chr> 
#>  1 X     block1
#>  2 X     block1
#>  3 Y     block1
#>  4 Y     block1
#>  5 Y     block1
#>  6 X     block2
#>  7 X     block2
#>  8 Y     block2
#>  9 Y     block2
#> 10 Y     block2
```

------------------------------------------------------------------------

## Method availability summary

| Method                                                                                         | multidesign | hyperdesign | multiframe | multiblock |
|------------------------------------------------------------------------------------------------|:-----------:|:-----------:|:----------:|:----------:|
| [`xdata()`](https://bbuchsbaum.github.io/multidesign/reference/xdata.md)                       |      ✓      |      ✓      |     ✓      |     \-     |
| [`design()`](https://bbuchsbaum.github.io/multidesign/reference/design.md)                     |      ✓      |      ✓      |     ✓      |     \-     |
| [`column_design()`](https://bbuchsbaum.github.io/multidesign/reference/column_design.md)       |      ✓      |      ✓      |     \-     |     \-     |
| [`subset()`](https://rdrr.io/r/base/subset.html)                                               |      ✓      |      ✓      |     ✓      |     \-     |
| [`split()`](https://rdrr.io/r/base/split.html)                                                 |      ✓      |     \-      |     ✓      |     \-     |
| [`select_variables()`](https://bbuchsbaum.github.io/multidesign/reference/select_variables.md) |      ✓      |      ✓      |     \-     |     \-     |
| [`summarize_by()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.md)         |      ✓      |     \-      |     ✓      |     \-     |
| [`fold_over()`](https://bbuchsbaum.github.io/multidesign/reference/fold_over.md)               |      ✓      |      ✓      |     ✓      |     \-     |
| [`cv_rows()`](https://bbuchsbaum.github.io/multidesign/reference/cv_rows.md)                   |      ✓      |      ✓      |     ✓      |     \-     |
| [`split_indices()`](https://bbuchsbaum.github.io/multidesign/reference/split_indices.md)       |      ✓      |     \-      |     ✓      |     \-     |
| [`block_indices()`](https://bbuchsbaum.github.io/multidesign/reference/block_indices.md)       |     \-      |      ✓      |     \-     |     ✓      |
| [`is_cstacked()`](https://bbuchsbaum.github.io/multidesign/reference/is_cstacked.md)           |     \-      |     \-      |     \-     |     ✓      |
| [`is_rstacked()`](https://bbuchsbaum.github.io/multidesign/reference/is_rstacked.md)           |     \-      |     \-      |     \-     |     ✓      |

------------------------------------------------------------------------

## Session info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.0       tibble_3.3.1      multidesign_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6         shape_1.4.6.1        xfun_0.57           
#>  [4] bslib_0.10.0         ggplot2_4.0.2        ggrepel_0.9.8       
#>  [7] lattice_0.22-9       vctrs_0.7.2          tools_4.5.3         
#> [10] generics_0.1.4       parallel_4.5.3       pkgconfig_2.0.3     
#> [13] multivarious_0.3.1   Matrix_1.7-4         RColorBrewer_1.1-3  
#> [16] S7_0.2.1             desc_1.4.3           assertthat_0.2.1    
#> [19] lifecycle_1.0.5      compiler_4.5.3       farver_2.1.2        
#> [22] GPArotation_2025.3-1 textshaping_1.0.5    codetools_0.2-20    
#> [25] htmltools_0.5.9      sass_0.4.10          yaml_2.3.12         
#> [28] deflist_0.2.0        glmnet_4.1-10        pillar_1.11.1       
#> [31] pkgdown_2.2.0        crayon_1.5.3         jquerylib_0.1.4     
#> [34] tidyr_1.3.2          cachem_1.1.0         iterators_1.0.14    
#> [37] foreach_1.5.2        parallelly_1.46.1    RSpectra_0.16-2     
#> [40] pls_2.9-0            svd_0.5.8            tidyselect_1.2.1    
#> [43] rsvd_1.0.5           digest_0.6.39        future_1.70.0       
#> [46] purrr_1.2.1          listenv_0.10.1       splines_4.5.3       
#> [49] fastmap_1.2.0        grid_4.5.3           cli_3.6.5           
#> [52] magrittr_2.0.4       utf8_1.2.6           survival_3.8-6      
#> [55] future.apply_1.20.2  withr_3.0.2          corpcor_1.6.10      
#> [58] scales_1.4.0         rmarkdown_2.31       globals_0.19.1      
#> [61] ragg_1.5.2           chk_0.10.0           memoise_2.0.1       
#> [64] evaluate_1.0.5       knitr_1.51           irlba_2.3.7         
#> [67] rlang_1.1.7          Rcpp_1.1.1           glue_1.8.0          
#> [70] geigen_2.3           jsonlite_2.0.0       R6_2.6.1            
#> [73] systemfonts_1.3.2    fs_2.0.1
```
