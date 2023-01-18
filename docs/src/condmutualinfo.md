# Conditional mutual information (CMI)

## API

The condition mutual information API is defined by

* [`ConditionalMutualInformation`](@ref),
* [`mutualinfo`](@ref),
* [`ConditionalMutualInformationEstimator`](@ref).

## Definitions

### Shannon CMI

```@docs
ConditionalMutualInformation
CMIShannon
CMIRenyiJizba
CMIRenyiPoczos
```

## Dedicated estimators

```@docs
condmutualinfo(::ConditionalMutualInformationEstimator, ::Any, ::Any, ::Any)
```

```@docs
ConditionalMutualInformationEstimator
FPVP
MesnerShalisi
PoczosSchneiderCMI
Rahimzamani
```

### [Table of dedicated CMI estimators](@id condmutualinfo_dedicated_estimators)

| Estimator                    | Principle         | [`CMIShannon`](@ref) | [`CMIRenyiPoczos`](@ref) |
| ---------------------------- | ----------------- | :------------------: | :----------------------: |
| [`FPVP`](@ref)               | Nearest neighbors |          ✓          |            x             |
| [`MesnerShalisi`](@ref)      | Nearest neighbors |          ✓          |            x             |
| [`Rahimzamani`](@ref)        | Nearest neighbors |          ✓          |            x             |
| [`PoczosSchneiderCMI`](@ref) | Nearest neighbors |          x           |            ✓            |

## Estimation through mutual information

```@docs
condmutualinfo(::MutualInformationEstimator, ::Any, ::Any, ::Any)
```

| Estimator                              |    Type    |     Principle     | [`CMIShannon`](@ref) |
| -------------------------------------- | :--------: | :---------------: | :------------------: |
| [`KraskovStögbauerGrassberger1`](@ref) | Continuous | Nearest neighbors |          ✓          |
| [`KraskovStögbauerGrassberger2`](@ref) | Continuous | Nearest neighbors |          ✓          |
| [`GaoKannanOhViswanath`](@ref)         |   Mixed    | Nearest neighbors |          ✓          |
| [`GaoOhViswanath`](@ref)               | Continuous | Nearest neighbors |          ✓          |

## Discrete CMI

```@docs
condmutualinfo(::ProbabilitiesEstimator, ::Any, ::Any, ::Any)
```

### [Table of discrete mutual information estimators](@id mutualinfo_overview)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`condmutualinfo`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           | [`CMIShannon`](@ref) | [`CMIRenyiSarbu`](@ref) |
| ---------------------------- | ------------------- | :------------------: | :---------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |          ✓          |           ✓            |
| [`ValueHistogram`](@ref)     | Binning (histogram) |          ✓          |           ✓            |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |          ✓          |           ✓            |
| [`Dispersion`](@ref)         | Dispersion patterns |          ✓          |           ✓            |

## Differential CMI

```@docs
condmutualinfo(::DifferentialEntropyEstimator, ::Any, ::Any, ::Any)
```

| Estimator                        | Principle         | Input data | [`CMIShannon`](@ref) | [`CMIRenyiSarbu`](@ref) |
| -------------------------------- | ----------------- | ---------- | :------------------: | :----------------: |
| [`Kraskov`](@ref)                | Nearest neighbors | `Dataset`  |          ✓          |         x          |
| [`Zhu`](@ref)                    | Nearest neighbors | `Dataset`  |          ✓          |         x          |
| [`ZhuSingh`](@ref)               | Nearest neighbors | `Dataset`  |          ✓          |         x          |
| [`Gao`](@ref)                    | Nearest neighbors | `Dataset`  |          ✓          |         x          |
| [`Goria`](@ref)                  | Nearest neighbors | `Dataset`  |          ✓          |         x          |
| [`Lord`](@ref)                   | Nearest neighbors | `Dataset`  |          ✓          |         x          |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors | `Dataset`  |          ✓          |         x          |
