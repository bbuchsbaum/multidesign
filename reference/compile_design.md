# Compile a design specification against an fmri frame

Compile a design specification against an fmri frame

## Usage

``` r
compile_design(frame, spec)
```

## Arguments

- frame:

  An \`fmri_frame\` or synchronized \`fmri_view\`.

- spec:

  A \`design_spec\`.

## Value

A \`compiled_design\` containing the dense fixed-effects model matrix,
component-aware term metadata, grouping data, and retained source spec.
