# Multivariate block formula special

\`mv()\` marks a named observation- or entity-aligned \`axis_block\` for
expansion by \[compile_design()\]. It is meaningful only inside a
formula.

## Usage

``` r
mv(block, components = NULL)
```

## Arguments

- block:

  Unquoted block name, optionally qualified by an entity name, such as
  \`motion\` or \`stimulus.visual_pca\`.

- components:

  Optional component positions or stable component IDs.

## Value

This function always errors when evaluated directly.
