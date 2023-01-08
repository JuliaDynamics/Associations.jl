# Conditional mutual information (CMI)

## API

```@docs
condmutualinfo
ConditionalMutualInformation
```

## Definitions

### Shannon CMI

```@docs
CMIShannon
CMIRenyiJizba
```

## Discrete

### [Table of discrete mutual information estimators](@id mutualinfo_overview)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`mutualinformation`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           |     Input data     | [`CMIShannon`](@ref) | [`CMIRenyiSarbu`](@ref) |
| ---------------------------- | ------------------- | :----------------: | :------------------: | :----------------: |
| [`CountOccurrences`](@ref)   | Frequencies         | `Vector`/`Dataset` |          ✅          |         ✅         |
| [`ValueHistogram`](@ref)     | Binning (histogram) | `Vector`/`Dataset` |          ✅          |         ✅         |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |      `Vector`      |          ✅          |         ✅         |
| [`Dispersion`](@ref)         | Dispersion patterns |      `Vector`      |          ✅          |         ✅         |

## Differential

### Estimators

The following are *differential mutual information* estimators. They always
override any definitions above with the concrete integral they estimate.

| Estimator                           | Principle         | [`CMIShannon`](@ref) | [`CMIRenyiSarbu`](@ref) |
| ----------------------------------- | ----------------- | :------------------: | :----------------: |
| [`FrenzelPompeVelmejkaPalus`](@ref) | Nearest neighbors |          ✅          |         x          |
| [`MesnerShalisi`](@ref)             | Nearest neighbors |          ✅          |         x          |
| [`PoczosSchneiderCMI`](@ref)        | Nearest neighbors |          x           |         ✅         |
| [`Rahimzamani`](@ref)               | Nearest neighbors |          ✅          |         x          |

Continuous/differential conditional mutual information may also be estimated using any of our
[`DifferentialEntropyEstimator`](@ref) that support multivariate input data.

| Estimator                        | Principle         | Input data | [`CMIShannon`](@ref) | [`CMIRenyiSarbu`](@ref) |
| -------------------------------- | ----------------- | ---------- | :------------------: | :----------------: |
| [`Kraskov`](@ref)                | Nearest neighbors | `Dataset`  |          ✅          |         x          |
| [`Zhu`](@ref)                    | Nearest neighbors | `Dataset`  |          ✅          |         x          |
| [`ZhuSingh`](@ref)               | Nearest neighbors | `Dataset`  |          ✅          |         x          |
| [`Gao`](@ref)                    | Nearest neighbors | `Dataset`  |          ✅          |         x          |
| [`Goria`](@ref)                  | Nearest neighbors | `Dataset`  |          ✅          |         x          |
| [`Lord`](@ref)                   | Nearest neighbors | `Dataset`  |          ✅          |         x          |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors | `Dataset`  |          ✅          |         x          |

```@docs
ConditionalMutualInformationEstimator
FrenzelPompeVelmejkaPalus
MesnerShalisi
PoczosSchneiderCMI
Rahimzamani
```
