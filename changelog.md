# Changelog

## 2.9

### Bug fixes

- Fixed bug in `transferentropy` function which yielded identical results in both directions for the bivariate case.
- Fixed bug that occurred when using `LocalPermutationTest` with `TEShannon` as the measure and a dedicated `TransferEntropyEstimator` (e.g. `Zhu1` or `Lindner`). This occurred because the `LocalPermutationTest` is, strictly speaking, a test using conditional mutual information as the measure. Therefore, naively applying a `TransferEntropy` measure such as `TEShannon` would error. This is fixed by performing a similar procedure where the source marginal is shuffled according to local neighborhoods in the conditional marginal. This is similar, but not identical to the CMI-based `LocalPermutationTest`, and adapts to the specific case of transfer entropy estimation using dedicated transfer entropy estimators instead of some lower-level estimator.
- Fixed bug in `Zhu1` transfer entropy estimator where when box volumes were extremely small, taking the logarithm of volume ratios resulted in `Inf` values. This was solved by simply ignoring these volumes.

## 2.8.0

Moved to DynamicalSystemsBase v3.0 (`trajectory` now returns both the data and the time
indices).

### Bugfixes

- Fixed bug in `GaussianMI` that occurred when the keyword `normalize` was set to `true`.

## 2.7.1

- Fixed an import warning.

## 2.7.0

- New association measure: `PMI` (part mutual information).

## 2.6.0

- New causal graph inference algorithm: `PC`.

## 2.5.0

- New independence test: `CorrTest`, based on (partial) correlations.

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
