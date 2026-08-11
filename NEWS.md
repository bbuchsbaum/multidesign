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
