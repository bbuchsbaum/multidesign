# multidesign 0.1.0.9000 (Development)

* Added `design_spec()` and the `mv()` formula special for frame-native design
  definitions.
* Added `compile_design()` with observation- and entity-block expansion,
  component-aware term metadata, factor-preserving dense parity, and random
  grouping metadata. Imaging assays remain owned by `fmridataset`.

## Entity correspondence and column-space contracts

- `hyperdesign()` can now declare an entity key with `id`, a column-space
  contract with `space = "common"` or `space = "block"`, or deliberately
  positional correspondence with `positional = TRUE`.
- `correspondence()` exposes the global entity universe, local-to-global row
  maps, block incidence, pairwise overlap counts, and connectivity. First-seen
  entity order remains the default; `sort_ids = TRUE` requests deterministic
  sorted order without mutating the source object.
- `entity_id()`, `column_space()`, and `has_correspondence()` expose the
  contract without requiring consumers to inspect attributes directly.
- Subsetting reconstructs the contract and recomputes correspondence.
  `select_variables()` revalidates common column spaces, while row folds retain
  ordinary local-row semantics.

## Explicit cell-observation masks

- `multidesign(..., cells =)` accepts an optional logical matrix matching the
  data dimensions. `cell_mask()` and `has_cell_mask()` query it. Missing values
  in `x` never create a mask implicitly, and an all-`FALSE` mask row remains a
  real design row.
- Row and column subsets, splits, folds, and variable selection propagate masks
  exactly. `bind_multidesign()` and `as_multidesign()` row-stack masks and treat
  unmasked inputs as all observed when another input has an explicit mask.
- Masked `summarize_by()` requires an explicit `aggregate` rule. Duplicate
  entity keys can likewise be collapsed explicitly in `hyperdesign()` with
  `aggregate = "mean"` or a scalar-returning function; conflicting non-key
  metadata remains an error.
- Preprocessing and dimensionality reduction reject partial masks until those
  operations can declare how observation masks transform.

## Explicit entity alignment

- `align_by_id()` joins common-space blocks onto the global entity universe and
  returns an `aligned_hyperdesign` containing a data array, a same-shaped cell
  mask, and a separate entity-by-block row-presence matrix.
- Absent rows and explicitly unobserved cells receive the requested display
  `fill`, but `observed` and `cells` remain authoritative. An observed `NA`, an
  unobserved cell, and an absent row therefore remain distinguishable.
- `as_multidesign()` remains a row-stack operation. `align_by_id()` is a data
  organization helper, not a generalized Procrustes solver; transforms,
  gauges, fitting, and certificates remain consumer responsibilities.

## Compatibility and maintenance

- Existing objects without correspondence or mask contracts retain their
  previous behavior.
- No new hard dependency was added. The unused `recipes` import was removed.
