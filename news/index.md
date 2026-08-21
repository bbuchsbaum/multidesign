# Changelog

## multidesign 0.1.0.9000 (Development)

- Added
  [`design_spec()`](https://bbuchsbaum.github.io/multidesign/reference/design_spec.md)
  and the
  [`mv()`](https://bbuchsbaum.github.io/multidesign/reference/mv.md)
  formula special for frame-native design definitions.
- Added
  [`compile_design()`](https://bbuchsbaum.github.io/multidesign/reference/compile_design.md)
  with observation- and entity-block expansion, component-aware term
  metadata, factor-preserving dense parity, and random grouping
  metadata. Imaging assays remain owned by `fmridataset`.

### Entity correspondence and column-space contracts

- [`hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.md)
  can now declare an entity key with `id`, a column-space contract with
  `space = "common"` or `space = "block"`, or deliberately positional
  correspondence with `positional = TRUE`.
- [`correspondence()`](https://bbuchsbaum.github.io/multidesign/reference/correspondence.md)
  exposes the global entity universe, local-to-global row maps, block
  incidence, pairwise overlap counts, and connectivity. First-seen
  entity order remains the default; `sort_ids = TRUE` requests
  deterministic sorted order without mutating the source object.
- [`entity_id()`](https://bbuchsbaum.github.io/multidesign/reference/entity_id.md),
  [`column_space()`](https://bbuchsbaum.github.io/multidesign/reference/column_space.md),
  and
  [`has_correspondence()`](https://bbuchsbaum.github.io/multidesign/reference/has_correspondence.md)
  expose the contract without requiring consumers to inspect attributes
  directly.
- Subsetting reconstructs the contract and recomputes correspondence.
  [`select_variables()`](https://bbuchsbaum.github.io/multidesign/reference/select_variables.md)
  revalidates common column spaces, while row folds retain ordinary
  local-row semantics.

### Explicit cell-observation masks

- `multidesign(..., cells =)` accepts an optional logical matrix
  matching the data dimensions.
  [`cell_mask()`](https://bbuchsbaum.github.io/multidesign/reference/cell_mask.md)
  and
  [`has_cell_mask()`](https://bbuchsbaum.github.io/multidesign/reference/cell_mask.md)
  query it. Missing values in `x` never create a mask implicitly, and an
  all-`FALSE` mask row remains a real design row.
- Row and column subsets, splits, folds, and variable selection
  propagate masks exactly.
  [`bind_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/bind_multidesign.md)
  and
  [`as_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/as_multidesign.md)
  row-stack masks and treat unmasked inputs as all observed when another
  input has an explicit mask.
- Masked
  [`summarize_by()`](https://bbuchsbaum.github.io/multidesign/reference/summarize_by.md)
  requires an explicit `aggregate` rule. Duplicate entity keys can
  likewise be collapsed explicitly in
  [`hyperdesign()`](https://bbuchsbaum.github.io/multidesign/reference/hyperdesign.md)
  with `aggregate = "mean"` or a scalar-returning function; conflicting
  non-key metadata remains an error.
- Preprocessing and dimensionality reduction reject partial masks until
  those operations can declare how observation masks transform.

### Explicit entity alignment

- [`align_by_id()`](https://bbuchsbaum.github.io/multidesign/reference/align_by_id.md)
  joins common-space blocks onto the global entity universe and returns
  an `aligned_hyperdesign` containing a data array, a same-shaped cell
  mask, and a separate entity-by-block row-presence matrix.
- Absent rows and explicitly unobserved cells receive the requested
  display `fill`, but `observed` and `cells` remain authoritative. An
  observed `NA`, an unobserved cell, and an absent row therefore remain
  distinguishable.
- [`as_multidesign()`](https://bbuchsbaum.github.io/multidesign/reference/as_multidesign.md)
  remains a row-stack operation.
  [`align_by_id()`](https://bbuchsbaum.github.io/multidesign/reference/align_by_id.md)
  is a data organization helper, not a generalized Procrustes solver;
  transforms, gauges, fitting, and certificates remain consumer
  responsibilities.

### Compatibility and maintenance

- Existing objects without correspondence or mask contracts retain their
  previous behavior.
- No new hard dependency was added. The unused `recipes` import was
  removed.
