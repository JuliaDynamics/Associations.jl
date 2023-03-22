# Changelog

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
