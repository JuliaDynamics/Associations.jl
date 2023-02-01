# Mutual information

## Mutual information API

The mutual information API is defined by

* [`MutualInformation`](@ref),
* [`mutualinfo`](@ref),
* [`MutualInformationEstimator`](@ref).

We provide a suite of estimators of various mutual information quantities. Many more
variants exist in the literature. Pull requests are welcome!

## Mutual information definitions

```@docs
MutualInformation
MIShannon
MITsallisFuruichi
MITsallisMartin
MIRenyiSarbu
MIRenyiJizba
```

## Dedicated estimators

```@docs
mutualinfo(est::MutualInformationEstimator, ::Any, ::Any) 
```

```@docs
MutualInformationEstimator
KraskovStögbauerGrassberger1
KraskovStögbauerGrassberger2
GaoKannanOhViswanath
GaoOhViswanath
```

### [Table of dedicated estimators](@id dedicated_estimators_mi)

| Estimator                              |    Type    |     Principle     | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiSarbu`](@ref) | [`MIRenyiJizba`](@ref) |
| -------------------------------------- | :--------: | :---------------: | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: | :--------------------: |
| [`KraskovStögbauerGrassberger1`](@ref) | Continuous | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`KraskovStögbauerGrassberger2`](@ref) | Continuous | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`GaoKannanOhViswanath`](@ref)         |   Mixed    | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`GaoOhViswanath`](@ref)               | Continuous | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |

## Discrete mutual information

```@docs
mutualinfo(::ProbabilitiesEstimator, ::Any, ::Any)
```

### [Table of discrete mutual information estimators](@id @id dedicated_probabilities_estimators_mi)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that can be used to compute discrete
[`mutualinformation`](@ref).

| Estimator                    | Principle           | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiJizba`](@ref) | [`MIRenyiSarbu`](@ref) |
| ---------------------------- | ------------------- | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: | :--------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |         ✓          |             ✓              |            ✓             |           ✓           |           x           |
| [`ValueHistogram`](@ref)     | Binning (histogram) |         ✓          |             ✓              |            ✓             |           ✓           |           x           |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |         ✓          |             ✓              |            ✓             |           ✓           |           x           |
| [`Dispersion`](@ref)         | Dispersion patterns |         ✓          |             ✓              |            ✓             |           ✓           |           x           |

### [Contingency matrix](@id contingency_matrix_mi)

```@docs
mutualinfo(::MutualInformation, ::ContingencyMatrix)
```

Discrete mutual information can be computed directly from its double-sum definition
by using the probabilities from a [`ContingencyMatrix`](@ref). This estimation
method works for    both numerical and categorical data, and the following
[`MutualInformation`](@ref)s are supported.

|                             | [`ContingencyMatrix`](@ref) |
| --------------------------- | :-------------------------: |
| [`MIShannon`](@ref)         |             ✓              |
| [`MITsallisFuruichi`](@ref) |             ✓              |
| [`MITsallisMartin`](@ref)   |             ✓              |
| [`MIRenyiSarbu`](@ref)      |             ✓              |
| [`MIRenyiJizba`](@ref)      |             ✓              |

## Differential/continuous mutual information

```@docs
mutualinfo(::DifferentialEntropyEstimator, ::Any, ::Any)
```

### [Table of differential mutual information estimators](@id dedicated_diffentropy_estimators_mi)

In addition to the dedicated differential mutual information estimators listed above,
continuous/differential mutual information may also be estimated using any of our
[`DifferentialEntropyEstimator`](@ref) that support multivariate input data.
When using these estimators, mutual information is computed as a sum
of entropy terms (with different dimensions), and no bias correction is applied.

| Estimator                        | Principle         | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiJizba`](@ref) | [`MIRenyiSurbu`](@ref) |
| -------------------------------- | ----------------- | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: | :--------------------: |
| [`Kraskov`](@ref)                | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Zhu`](@ref)                    | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`ZhuSingh`](@ref)               | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Gao`](@ref)                    | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Goria`](@ref)                  | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Lord`](@ref)                   | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
