# Getting Started with multidesign

``` r
if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("albersdown", quietly = TRUE)) ggplot2::theme_set(albersdown::theme_albers(family = params$family, preset = params$preset))
library(multidesign)
library(tibble)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Why multidesign?

In experimental research—particularly neuroimaging, psychophysics, and
behavioral studies—data often has a specific structure that standard R
data frames don’t capture well:

- **Observations** (trials, time points, scans) arranged as rows in a
  matrix
- **Variables** (brain regions, sensors, features) arranged as columns
- **Design metadata** (condition, subject, run, block) describing each
  observation
- **Variable metadata** (ROI name, hemisphere, sensor type) describing
  each column

Without a unified structure, you end up with separate objects for data
and metadata that can easily get out of sync, manual bookkeeping for
cross-validation splits, and difficulty handling multi-subject or
multi-session experiments.

The `multidesign` package solves these problems by providing data
structures that:

1.  **Keep data and metadata together** — No more matching row indices
    between separate objects
2.  **Support design-aware operations** — Split, subset, and summarize
    by experimental factors
3.  **Handle multi-block data naturally** — Multiple subjects, sessions,
    or modalities
4.  **Provide built-in cross-validation** — Proper stratification
    respecting your experimental design

## Quick reference: When to use each class

| Class         | Use when…                                                                                             |
|---------------|-------------------------------------------------------------------------------------------------------|
| `multidesign` | You have a single data matrix with row-wise design information and optional column metadata           |
| `hyperdesign` | You have multiple related matrices (e.g., different subjects or sessions) that share design structure |
| `multiframe`  | Your data is expensive to load and you want lazy evaluation (load on demand)                          |
| `multiblock`  | You have stacked matrices sharing a row or column dimension (low-level operations)                    |

Most users will start with `multidesign` and move to `hyperdesign` when
working with multi-subject data.

## A complete workflow

Let’s walk through a realistic analysis workflow using simulated
neuroimaging-style data.

### Step 1: Create your data

Imagine a simple experiment: 2 conditions (face vs. house), 10 subjects,
20 trials per condition per subject, measuring activity in 50 brain
regions.

``` r
# Simulate data for one subject
n_trials <- 40
n_regions <- 50

# Create the data matrix (trials x brain regions)
X <- matrix(rnorm(n_trials * n_regions), n_trials, n_regions)

# Create the design information
# Each run has 10 trials: 5 face + 5 house
design_df <- tibble(
 condition = rep(rep(c("face", "house"), each = 5), 4),
 trial = 1:40,
 run = rep(1:4, each = 10)
)

# Create column metadata describing the brain regions
column_info <- tibble(
 region = paste0("ROI_", 1:n_regions),
 hemisphere = rep(c("left", "right"), 25),
 network = rep(c("visual", "frontal", "parietal", "temporal", "default"), 10)
)
```

### Step 2: Build a multidesign

Combine your data matrix, design information, and column metadata into a
single object:

``` r
mds <- multidesign(X, design_df, column_info)
mds
#> 
#> === Multidesign Object ===
#> 
#> Data Matrix: 
#>    40 observations x 50 variables 
#> 
#> Design Variables: 
#>   * condition: 2 levels (face, house)
#>   * trial: 40 levels (1, 2, 3...39, 40)
#>   * run: 4 levels (1, 2, 3, 4)
#> 
#> Column Metadata:
#>   * region: 50 levels (ROI_1, ROI_2, ROI_3...ROI_49, ROI_50)
#>   * hemisphere: 2 levels (left, right)
#>   * network: 5 levels (visual, frontal, parietal, temporal, default)
#> 
#> =======================
#> 
```

The `multidesign` object keeps everything synchronized. You can extract
components with accessor functions:

``` r
# Get the data matrix
head(xdata(mds)[, 1:5])
#>            [,1]       [,2]        [,3]         [,4]       [,5]
#> [1,]  1.3709584  0.2059986  1.51270701 -1.493625067 -0.1755259
#> [2,] -0.5646982 -0.3610573  0.25792144 -1.470435741 -1.0717824
#> [3,]  0.3631284  0.7581632  0.08844023  0.124702386  0.1632069
#> [4,]  0.6328626 -0.7267048 -0.12089654 -0.996639135 -0.3627384
#> [5,]  0.4042683 -1.3682810 -1.19432890 -0.001822614  0.5900135
#> [6,] -0.1061245  0.4328180  0.61199690 -0.428258881  1.4324219

# Get the design (without internal .index column)
head(design(mds))
#> # A tibble: 6 × 3
#>   condition trial   run
#>   <chr>     <int> <int>
#> 1 face          1     1
#> 2 face          2     1
#> 3 face          3     1
#> 4 face          4     1
#> 5 face          5     1
#> 6 house         6     1

# Get column metadata
head(column_design(mds))
#> # A tibble: 6 × 3
#>   region hemisphere network 
#>   <chr>  <chr>      <chr>   
#> 1 ROI_1  left       visual  
#> 2 ROI_2  right      frontal 
#> 3 ROI_3  left       parietal
#> 4 ROI_4  right      temporal
#> 5 ROI_5  left       default 
#> 6 ROI_6  right      visual
```

### Step 3: Explore and subset your data

Filter observations by design variables:

``` r
# Keep only face trials
mds_face <- subset(mds, condition == "face")
mds_face
#> 
#> === Multidesign Object ===
#> 
#> Data Matrix: 
#>    20 observations x 50 variables 
#> 
#> Design Variables: 
#>   * condition: 1 levels (face)
#>   * trial: 20 levels (1, 2, 3...34, 35)
#>   * run: 4 levels (1, 2, 3, 4)
#> 
#> Column Metadata:
#>   * region: 50 levels (ROI_1, ROI_2, ROI_3...ROI_49, ROI_50)
#>   * hemisphere: 2 levels (left, right)
#>   * network: 5 levels (visual, frontal, parietal, temporal, default)
#> 
#> =======================
#> 

# Keep only run 1
mds_run1 <- subset(mds, run == 1)
nrow(xdata(mds_run1))
#> [1] 10
```

Select variables (columns) by their metadata:

``` r
# Keep only visual network regions
mds_visual <- select_variables(mds, network == "visual")
ncol(xdata(mds_visual))
#> [1] 10
column_design(mds_visual)
#> # A tibble: 10 × 3
#>    region hemisphere network
#>    <chr>  <chr>      <chr>  
#>  1 ROI_1  left       visual 
#>  2 ROI_6  right      visual 
#>  3 ROI_11 left       visual 
#>  4 ROI_16 right      visual 
#>  5 ROI_21 left       visual 
#>  6 ROI_26 right      visual 
#>  7 ROI_31 left       visual 
#>  8 ROI_36 right      visual 
#>  9 ROI_41 left       visual 
#> 10 ROI_46 right      visual

# Keep only left hemisphere visual regions
mds_left_visual <- select_variables(mds, network == "visual" & hemisphere == "left")
ncol(xdata(mds_left_visual))
#> [1] 5
```

### Step 4: Summarize by design factors

Compute summary statistics grouped by design variables:

``` r
# Mean activity per condition
mds_means <- summarize_by(mds, condition)
mds_means
#> 
#> === Multidesign Object ===
#> 
#> Data Matrix: 
#>    2 observations x 50 variables 
#> 
#> Design Variables: 
#>   * condition: 2 levels (face, house)
#> 
#> Column Metadata:
#>   * region: 50 levels (ROI_1, ROI_2, ROI_3...ROI_49, ROI_50)
#>   * hemisphere: 2 levels (left, right)
#>   * network: 5 levels (visual, frontal, parietal, temporal, default)
#> 
#> =======================
#> 

# The result is a new multidesign with one row per group
dim(xdata(mds_means))
#> [1]  2 50
design(mds_means)
#> # A tibble: 2 × 1
#>   condition
#>   <chr>    
#> 1 face     
#> 2 house
```

### Step 5: Split for detailed analysis

Split into separate multidesign objects:

``` r
# Split by condition
by_condition <- split(mds, condition)
length(by_condition)
#> [1] 2
names(by_condition) <- c("face", "house")

# Each element is a complete multidesign
by_condition$face
#> 
#> === Multidesign Object ===
#> 
#> Data Matrix: 
#>    20 observations x 50 variables 
#> 
#> Design Variables: 
#>   * condition: 1 levels (face)
#>   * trial: 20 levels (1, 2, 3...34, 35)
#>   * run: 4 levels (1, 2, 3, 4)
#> 
#> Column Metadata:
#>   * region: 50 levels (ROI_1, ROI_2, ROI_3...ROI_49, ROI_50)
#>   * hemisphere: 2 levels (left, right)
#>   * network: 5 levels (visual, frontal, parietal, temporal, default)
#> 
#> =======================
#> 
```

### Step 6: Cross-validation

Create cross-validation folds that respect your experimental design:

``` r
# Leave-one-run-out cross-validation
folds <- fold_over(mds, run)
folds
#> 
#> === Cross-Validation Folds ===
#> 
#>  Number of Folds: 4 
#> 
#> Fold 1: 
#>   * Analysis Set: 30 observations
#>   * Assessment Set: 10 observations
#> 
#> Fold 2: 
#>   * Analysis Set: 30 observations
#>   * Assessment Set: 10 observations
#> 
#> Fold 3: 
#>   * Analysis Set: 30 observations
#>   * Assessment Set: 10 observations
#> 
#> Fold 4: 
#>   * Analysis Set: 30 observations
#>   * Assessment Set: 10 observations
#> 
#> =============================
#> 

# Each fold contains analysis (training) and assessment (test) sets
fold1 <- folds[[1]]
nrow(xdata(fold1$analysis))
#> [1] 30
nrow(xdata(fold1$assessment))
#> [1] 10
```

If downstream code needs stable source-row identities, ask fold
construction to carry them through:

``` r
folds_with_ids <- fold_over(mds, run, preserve_row_ids = TRUE)
folds_with_ids[[1]]$assessment$design$.orig_index
#>  [1]  1  2  3  4  5  6  7  8  9 10
```

### Step 7: Run a cross-validation pipeline

Use
[`cross_validate()`](https://bbuchsbaum.github.io/multidesign/reference/cross_validate.md)
to execute a complete CV workflow:

``` r
# Define a simple "model": compute condition means on training data
fit_fn <- function(analysis) {
 means_by_cond <- summarize_by(analysis, condition)
 list(
   face_mean = xdata(means_by_cond)[1, ],
   house_mean = xdata(means_by_cond)[2, ]
 )
}

# Define scoring: correlation between predicted and actual means
score_fn <- function(model, assessment) {
 test_means <- summarize_by(assessment, condition)

 # Compare model predictions to test data
 face_cor <- cor(model$face_mean, xdata(test_means)[1, ])
 house_cor <- cor(model$house_mean, xdata(test_means)[2, ])

 c(face_r = face_cor, house_r = house_cor, mean_r = (face_cor + house_cor) / 2)
}

# Run cross-validation
cv_result <- cross_validate(folds, fit_fn, score_fn)
cv_result
#> 
#> === Cross-Validation Result ===
#> 
#> Folds: 4 
#> 
#> Metrics:
#>    face_r : mean = 0.0176 , sd = 0.1051 
#>    house_r : mean = 0.0332 , sd = 0.0865 
#>    mean_r : mean = 0.0254 , sd = 0.0626 
#> 
#> ===============================
#> 
```

Examine detailed results:

``` r
# Get summary statistics
summary(cv_result)
#> # A tibble: 3 × 6
#>   metric    mean     sd     min    max  median
#>   <chr>    <dbl>  <dbl>   <dbl>  <dbl>   <dbl>
#> 1 face_r  0.0176 0.105  -0.100  0.153  0.00883
#> 2 house_r 0.0332 0.0865 -0.0891 0.113  0.0546 
#> 3 mean_r  0.0254 0.0626 -0.0497 0.0982 0.0266

# Access raw scores
cv_result$scores
#> # A tibble: 4 × 4
#>    face_r house_r   mean_r .fold
#>     <dbl>   <dbl>    <dbl> <int>
#> 1 -0.100   0.113   0.00624     1
#> 2 -0.0103 -0.0891 -0.0497      2
#> 3  0.153   0.0432  0.0982      3
#> 4  0.0280  0.0659  0.0470      4
```

## Working with multiple subjects: hyperdesign

Real experiments typically have multiple subjects. A `hyperdesign`
manages multiple `multidesign` objects that share common design
structure.

``` r
# Create data for 3 subjects
subjects <- lapply(1:3, function(subj) {
 X_subj <- matrix(rnorm(n_trials * n_regions), n_trials, n_regions)

 # Add a subject-specific signal
 if (subj == 1) X_subj[1:20, 1:10] <- X_subj[1:20, 1:10] + 0.5
 if (subj == 2) X_subj[1:20, 1:10] <- X_subj[1:20, 1:10] + 0.3
 if (subj == 3) X_subj[1:20, 1:10] <- X_subj[1:20, 1:10] + 0.7

 multidesign(X_subj, design_df, column_info)
})

# Combine into a hyperdesign
hd <- hyperdesign(subjects, block_names = c("subj1", "subj2", "subj3"))
hd
#> 
#> === Hyperdesign Object ===
#> 
#> Number of blocks:  3 
#> 
#> +- Block  1  (subj1)  -----------------
#> |  Dimensions: 40 x 50 
#> |  Design Variables: condition, trial, run 
#> |  Design Structure: 
#> |   * condition: 2 levels (face, house)
#> |   * trial: 40 levels (1, 2, 3...39, 40)
#> |   * run: 4 levels (1, 2, 3, 4)
#> |  Column Design: Present
#> |   Variables:  region, hemisphere, network 
#> 
#> +- Block  2  (subj2)  -----------------
#> |  Dimensions: 40 x 50 
#> |  Design Variables: condition, trial, run 
#> |  Design Structure: 
#> |   * condition: 2 levels (face, house)
#> |   * trial: 40 levels (1, 2, 3...39, 40)
#> |   * run: 4 levels (1, 2, 3, 4)
#> |  Column Design: Present
#> |   Variables:  region, hemisphere, network 
#> 
#> +- Block  3  (subj3)  -----------------
#> |  Dimensions: 40 x 50 
#> |  Design Variables: condition, trial, run 
#> |  Design Structure: 
#> |   * condition: 2 levels (face, house)
#> |   * trial: 40 levels (1, 2, 3...39, 40)
#> |   * run: 4 levels (1, 2, 3, 4)
#> |  Column Design: Present
#> |   Variables:  region, hemisphere, network 
#> 
#> =======================
#> 
```

### Leave-one-subject-out cross-validation

``` r
# Without specifying variables, fold_over creates leave-one-block-out folds
loso_folds <- fold_over(hd)
length(loso_folds)
#> [1] 3

# Each fold: train on 2 subjects, test on 1
fold1 <- loso_folds[[1]]
length(fold1$analysis)
#> [1] 2
fold1$assessment
#> 
#> === Multidesign Object ===
#> 
#> Data Matrix: 
#>    40 observations x 50 variables 
#> 
#> Design Variables: 
#>   * condition: 2 levels (face, house)
#>   * trial: 40 levels (1, 2, 3...39, 40)
#>   * run: 4 levels (1, 2, 3, 4)
#> 
#> Column Metadata:
#>   * region: 50 levels (ROI_1, ROI_2, ROI_3...ROI_49, ROI_50)
#>   * hemisphere: 2 levels (left, right)
#>   * network: 5 levels (visual, frontal, parietal, temporal, default)
#> 
#> =======================
#> 
```

### Within-subject cross-validation

``` r
# Fold by run within each subject
run_folds <- fold_over(hd, run)
length(run_folds)
#> [1] 12
```

### Explicit row-based cross-validation

When you already know the held-out rows, use
[`cv_rows()`](https://bbuchsbaum.github.io/multidesign/reference/cv_rows.md)
to build a `foldlist` directly:

``` r
row_folds <- cv_rows(
  hd,
  rows = list(
    list(subj1 = 1:2, subj2 = 1:2),
    list(subj1 = 3:4, subj2 = 3:4)
  )
)

row_folds[[1]]$assessment
#> 
#> === Hyperdesign Object ===
#> 
#> Number of blocks:  2 
#> 
#> +- Block  1  (subj1)  -----------------
#> |  Dimensions: 2 x 50 
#> |  Design Variables: condition, trial, run 
#> |  Design Structure: 
#> |   * condition: 1 levels (face)
#> |   * trial: 2 levels (1, 2)
#> |   * run: 1 levels (1)
#> |  Column Design: Present
#> |   Variables:  region, hemisphere, network 
#> 
#> +- Block  2  (subj2)  -----------------
#> |  Dimensions: 2 x 50 
#> |  Design Variables: condition, trial, run 
#> |  Design Structure: 
#> |   * condition: 1 levels (face)
#> |   * trial: 2 levels (1, 2)
#> |   * run: 1 levels (1)
#> |  Column Design: Present
#> |   Variables:  region, hemisphere, network 
#> 
#> =======================
#> 
```

Use `preserve_row_ids = TRUE` here as well if external correspondence
tables or annotations are keyed by original rows:

``` r
row_folds_with_ids <- cv_rows(
  hd,
  rows = list(list(subj1 = 1:2, subj2 = 1:2)),
  preserve_row_ids = TRUE
)

row_folds_with_ids[[1]]$assessment[[1]]$design$.orig_index
#> [1] 1 2
row_folds_with_ids[[1]]$held_out$row_ids
#> $subj1
#> [1] 1 2
#> 
#> $subj2
#> [1] 1 2
```

### Collapse to a single multidesign

If you need to analyze all subjects together:

``` r
# Collapse hyperdesign to one multidesign
mds_all <- as_multidesign(hd, .id = "subject")
mds_all
#> 
#> === Multidesign Object ===
#> 
#> Data Matrix: 
#>    120 observations x 50 variables 
#> 
#> Design Variables: 
#>   * condition: 2 levels (face, house)
#>   * trial: 40 levels (1, 2, 3...39, 40)
#>   * run: 4 levels (1, 2, 3, 4)
#>   * subject: 3 levels (subj1, subj2, subj3)
#> 
#> Column Metadata:
#>   * region: 50 levels (ROI_1, ROI_2, ROI_3...ROI_49, ROI_50)
#>   * hemisphere: 2 levels (left, right)
#>   * network: 5 levels (visual, frontal, parietal, temporal, default)
#> 
#> =======================
#> 

# The .id column tracks which block each observation came from
table(design(mds_all)$subject)
#> 
#> subj1 subj2 subj3 
#>    40    40    40
```

## Next steps

- See
  [`vignette("class_overview")`](https://bbuchsbaum.github.io/multidesign/articles/class_overview.md)
  for detailed documentation of each class
- Explore the function reference for methods like
  [`bind_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/bind_multidesign.md),
  [`df_to_hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/df_to_hyperdesign.md),
  and `reduce()`
- Check out the package website for additional examples

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
#> [58] scales_1.4.0         rmarkdown_2.31       albersdown_1.0.0    
#> [61] globals_0.19.1       ragg_1.5.2           chk_0.10.0          
#> [64] memoise_2.0.1        evaluate_1.0.5       knitr_1.51          
#> [67] irlba_2.3.7          rlang_1.1.7          Rcpp_1.1.1          
#> [70] glue_1.8.0           geigen_2.3           jsonlite_2.0.0      
#> [73] R6_2.6.1             systemfonts_1.3.2    fs_2.0.1
```
