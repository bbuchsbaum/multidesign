# multidesign 0.1.0.9000 (Development)

* Added `design_spec()` and the `mv()` formula special for frame-native design
  definitions.
* Added `compile_design()` with observation- and entity-block expansion,
  component-aware term metadata, factor-preserving dense parity, and random
  grouping metadata. Imaging assays remain owned by `fmridataset`.
* Completed relation-aware formula compilation: resolved entity scalars and
  lazily lifted entity blocks now share the canonical observation view,
  `mv()` syntax and component selectors are validated strictly, interaction
  attribution uses exact generated variables, and random formulas fail early
  on unsupported blocks or missing variables.
* Added reusable compiled-design blueprints. `apply_design()` now preserves
  training factor levels, explicit contrast matrices, formula transformation
  parameters, model columns, and stable `mv()` component IDs. `design_rows()`
  records fixed, random, and multivariate missingness under explicit fail or
  omit policies.
* Expanded compiled-design metadata into normalized, machine-readable tables.
  `term_data()` gives stable fixed-term grouping, `component_data()` records
  every model-column-to-block-component link, and `grouping_term_data()` plus
  `random_effect_data()` describe multiple `|` and `||` terms without requiring
  downstream formula parsing. Random-effect transformations now participate in
  the compiled missing-value policy.
* Added semantic `design_input_digest()` and `design_digest()` fingerprints,
  bounded mutation-isolated runtime caches for `compile_design()`, and
  `compile_design_folds()` for leakage-safe frame folds. Cache keys track only
  design dependencies and referenced block/component selectors, never imaging
  assays or feature layout; assessment designs always reuse an analysis-fitted
  blueprint.
* Added randomized dense-reference, permutation, backend, identifier, cache,
  and selector stress tests. Overlapping `mv()` terms now reuse one generated
  column per stable block/component identity, and generated names are fixed by
  the block's complete component domain so syntactically hostile IDs remain
  unambiguous across different selections.
