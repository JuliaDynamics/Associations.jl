# Conditional entropy

## Conditional entropy API

The mutual information API is defined by

* [`ConditionalEntropy`](@ref),
* [`entropy_conditional`](@ref),

We provide a suite of estimators of various mutual information quantities. Many more
variants exist in the literature. Pull requests are welcome!

## Conditional entropy definitions

```@docs
ConditionalEntropy
CEShannon
CETsallisFuruichi
CETsallisAbe
```

More variants exist in the literature. Pull requests are welcome!

## Discrete conditional entropy

```@docs
entropy_conditional(::ConditionalEntropy, ::ContingencyMatrix)
```

### [Contingency matrix](@id contingency_matrix_ce)

Discrete conditional entropy can be computed directly from its sum-definition
by using the probabilities from a [`ContingencyMatrix`](@ref). This estimation
method works for  both numerical and categorical data, and the following
[`ConditionalEntropy`](@ref) definitions are supported.

|                             | [`ContingencyMatrix`](@ref) |
| --------------------------- | :-------------------------: |
| [`CEShannon`](@ref)         |             ✓              |
| [`CETsallisFuruichi`](@ref) |             ✓              |
| [`CETsallisAbe`](@ref)      |             ✓              |

### [Table of discrete conditional entropy estimators](@id probabilities_estimators_ce)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`entropy_conditional`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           | [`CEShannon`](@ref) | [`CETsallisAbe`](@ref) | [`CETsallisFuruichi`](@ref) |
| ---------------------------- | ------------------- | :-----------------: | :--------------------: | :-------------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |         ✓          |           ✓           |              x              |
| [`ValueHistogram`](@ref)     | Binning (histogram) |         ✓          |           ✓           |              x              |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |         ✓          |           ✓           |              x              |
| [`Dispersion`](@ref)         | Dispersion patterns |         ✓          |           ✓           |              x              |

## Differential/continuous conditional entropy

### [Table of differential conditional entropy estimators](@id diffentropy_estimators_ce)

Continuous/differential mutual information may be estimated using any of our
[`DifferentialEntropyEstimator`](@ref)s that support multivariate input data.

| Estimator                        | Principle         | [`CEShannon`](@ref) | [`CETsallisAbe`](@ref) | [`CETsallisFuruichi`](@ref) |
| -------------------------------- | ----------------- | :-----------------: | :--------------------: | :-------------------------: |
| [`Kraskov`](@ref)                | Nearest neighbors |         ✓          |           x           |              x              |
| [`Zhu`](@ref)                    | Nearest neighbors |         ✓          |           x           |              x              |
| [`ZhuSingh`](@ref)               | Nearest neighbors |         ✓          |           x           |              x              |
| [`Gao`](@ref)                    | Nearest neighbors |         ✓          |           x           |              x              |
| [`Goria`](@ref)                  | Nearest neighbors |         ✓          |           x           |              x              |
| [`Lord`](@ref)                   | Nearest neighbors |         ✓          |           x           |              x              |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors |         ✓          |           x           |              x              |
