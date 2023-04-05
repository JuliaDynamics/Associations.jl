# Changelog

## 2.4.0

- Added partial distance correlation measure. To compute it, simply provide a
    third input argument to `distance_correlation`.
- `DistanceCorrelation` is now compatible with both `SurrogateTest` and
    `LocalPermutationTest` in its conditional form.

## 2.3.1

- The `MesnerShalisi` estimator is now deprecated and renamed to `MesnerShalizi` (with
    correct spelling).

## 2.3.0

- Significant speed-ups for `OCE` by sorting on maximal measure, thus avoiding
    unnecessary significance tests.
- Default parameters for `OCE` default lag parameter have changed. Now, `τmax = 1`, since
    that is the only case considered in the original paper. We also use the
    `MesnerShalisi` CMI estimator for the conditional step, because in contrast to
    the `FPVP` estimator, it has been shown to be consistent.
- Source code for `OCE` has been drastically simplified by merging the pairwise
    and conditional parent finding steps.
- `OCE` result can now be converted to a `SimpleDiGraph` from Graphs.jl.

## 2.2.1

- `infer_graph` now accepts the `verbose` keyword.
- Fixed a bug in backwards elimination step of `OCE` algorithm that was caused due
    to an undefined variable.

## 2.2

- Added `MCR` and `RMCD` recurrence based association measures, along with 
    the corresponding `mcr` and `rmcd` methods.

## 2.1

- Bugfix for `OCE` for certain conditional variable cases.
- Improve docstring for `OCE`.
- Updated readme.
- Fixed bug related to `DifferentialEntropyEstimator` "unit" conversion.

## 2.0

The entire package has been completely redesigned. Previous methods are deprecated,
and will disappear completely in a v2.X release.
For a full overview of new functionality, see the online documentation.

## 1.4.0

### New methods

- Added [`conditional_mutualinfo`](@ref).

## 1.3.0

### Example systems

- Added continous example system [`repressilator`](@ref).
