# multidesign

`multidesign` is an R package for multivariate experimental data where
the matrix alone is not enough. It keeps observation-level design
metadata and variable-level metadata attached to the data, so
subsetting, grouping, splitting, and cross-validation stay aligned.

This is especially useful for workflows in neuroimaging, behavioral
science, psychophysics, and other settings where you need to manage:

- observations such as trials, scans, or time points
- variables such as ROIs, sensors, or engineered features
- row-wise design factors such as condition, run, session, or subject
- multi-block datasets spanning subjects, sessions, or modalities

## Why use `multidesign`?

- Keep the data matrix, row design, and optional column metadata in one
  object.
- Filter observations with design-aware
  [`subset()`](https://rdrr.io/r/base/subset.html) semantics.
- Select variables using column metadata via
  [`select_variables()`](https://bbuchsbaum.github.io/multidesign/reference/select_variables.md).
- Split and summarize by experimental factors with
  [`split()`](https://rdrr.io/r/base/split.html) and
  [`summarize_by()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.md).
- Build cross-validation folds that respect your design with
  [`fold_over()`](https://bbuchsbaum.github.io/multidesign/reference/fold_over.md)
  and
  [`cv_rows()`](https://bbuchsbaum.github.io/multidesign/reference/cv_rows.md).
- Scale from a single matrix to multi-subject or lazy-loading workflows
  with `hyperdesign` and `multiframe`.

## Installation

`multidesign` is not on CRAN yet. Install the development version from
GitHub:

``` r
install.packages("remotes")
remotes::install_github("bbuchsbaum/multidesign")
```

## Core objects

| Class         | Use it when…                                                                            |
|---------------|-----------------------------------------------------------------------------------------|
| `multidesign` | You have one data matrix plus row-wise design information and optional column metadata. |
| `hyperdesign` | You have multiple related `multidesign` objects, such as one per subject or session.    |
| `multiframe`  | Your observations are expensive to materialize and you want lazy evaluation.            |
| `multiblock`  | You need lower-level stacked matrix operations without the full design abstraction.     |

Most workflows start with `multidesign` and move to `hyperdesign` when
the data naturally breaks into subject-, session-, or modality-level
blocks.

## Example

``` r
library(multidesign)

set.seed(42)

# Simulate a trial-by-feature matrix
X <- matrix(rnorm(40 * 8), nrow = 40, ncol = 8)

# Design metadata for each observation
design_df <- tibble::tibble(
  condition = rep(c("face", "house"), each = 20),
  run = rep(1:4, each = 10),
  trial = 1:40
)

# Metadata for each variable / column
column_df <- tibble::tibble(
  roi = paste0("ROI_", 1:8),
  network = rep(c("visual", "frontal"), each = 4)
)

# Build a multidesign object
md <- multidesign(X, design_df, column_df)

# Row-wise and column-wise access stay synchronized
face_trials <- subset(md, condition == "face")
visual_rois <- select_variables(md, network == "visual")

# Summarize by an experimental factor
condition_means <- summarize_by(md, condition)

# Create design-aware cross-validation folds
folds <- fold_over(md, run)

# Run a simple CV pipeline
cv <- cross_validate(
  folds,
  fit_fn = function(analysis) {
    list(mean = colMeans(xdata(analysis)))
  },
  score_fn = function(model, assessment) {
    mat <- xdata(assessment)
    pred <- matrix(
      model$mean,
      nrow = nrow(mat),
      ncol = ncol(mat),
      byrow = TRUE
    )
    c(mse = mean((mat - pred)^2))
  }
)

summary(cv)
```

## Typical workflow

1.  Create a `multidesign` from a matrix and a design table.
2.  Subset rows by design variables and columns by variable metadata.
3.  Summarize or split by condition, run, subject, or any other design
    factor.
4.  Generate folds with
    [`fold_over()`](https://bbuchsbaum.github.io/multidesign/reference/fold_over.md)
    or
    [`cv_rows()`](https://bbuchsbaum.github.io/multidesign/reference/cv_rows.md).
5.  Pass those folds into
    [`cross_validate()`](https://bbuchsbaum.github.io/multidesign/reference/cross_validate.md)
    with your model-fitting and scoring functions.
6.  Move to `hyperdesign` when the analysis spans multiple related
    blocks.

## Learn more

- [Package website](https://bbuchsbaum.github.io/multidesign/)
- [Getting started
  vignette](https://bbuchsbaum.github.io/multidesign/vignettes/Introduction.Rmd)
- [Class overview
  vignette](https://bbuchsbaum.github.io/multidesign/vignettes/class_overview.Rmd)
- [Issue tracker](https://github.com/bbuchsbaum/multidesign/issues)

## Development status

The package already includes tests for the core data structures and
cross-validation pipeline, but the API is still evolving. Expect the
README and vignettes to be the best entry points while the package
matures.

## Albers theme

This package uses the albersdown theme. Existing vignette theme hooks
are replaced so `albers.css` and local `albers.js` render consistently
on CRAN and GitHub Pages. The defaults are configured via
`params$family` and `params$preset` (family = “ochre”, preset =
“homage”). The pkgdown site uses `template: { package: albersdown }`.
