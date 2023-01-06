# Mutual information

The function that estimates mutual information from data is [`mutualinfo`](@ref).
It can estimate different types of mutual informations from the scientific literature,
and each is represented here as a subtype of [`MutualInformation`](@ref). Because a
mutual information can be formulated in many different ways, each mutual information
type can be estimated according to multiple definition, which are represented by
subtypes of [`MutualInformationDefinition`](@ref)).

To see which estimators are compatible with the various definitions, see the
[overview table](@ref mutualinfo_overview) below.

## API

```@docs
mutualinfo
MutualInformation
MutualInformationDefinition
```

## Types of mutual information

### Shannon mutual information

```@docs
MIShannon
```

### Tsallis mutual information

```@docs
MITsallisFuruichi
MITsallisMartin
```

### Rényi mutual information

```@docs
MIRenyiSarbu
```

More variants of the discrete Renyi mutual information are possible.
Pull requests are welcome!

## Discrete mutual information

### [Table of discrete mutual information estimators](@id mutualinfo_overview)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`mutualinformation`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           |     Input data     | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiSurbu`](@ref) |
| ---------------------------- | ------------------- | :----------------: | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         | `Vector`/`Dataset` |         ✅          |             ✅              |            ✅             |           ✅           |
| [`ValueHistogram`](@ref)     | Binning (histogram) | `Vector`/`Dataset` |         ✅          |             ✅              |            ✅             |           ✅           |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |      `Vector`      |         ✅          |             ✅              |            ✅             |           ✅           |
| [`Dispersion`](@ref)         | Dispersion patterns |      `Vector`      |         ✅          |             ✅              |            ✅             |           ✅           |

## Differential mutual information

### Estimators

The following are *differential mutual information* estimators.

| Estimator                              | Principle         | [`MIShannon`](@ref) |
| -------------------------------------- | ----------------- | :-----------------: |
| [`KraskovStögbauerGrassberger1`](@ref) | Nearest neighbors |         ✅          |
| [`KraskovStögbauerGrassberger2`](@ref) | Nearest neighbors |         ✅          |
| [`GaoKannanOhViswanath`](@ref)         | Nearest neighbors |         ✅          |
| [`GaoOhViswanath`](@ref)               | Nearest neighbors |         ✅          |

Continuous/differential mutual information may also be estimated using any of our
[`DifferentialEntropyEstimator`](@ref) that support multivariate input data.

| Estimator                        | Principle         | Input data | [`MIShannon`](@ref) |
| -------------------------------- | ----------------- | ---------- | :-----------------: |
| [`Kraskov`](@ref)                | Nearest neighbors | `Dataset`  |         ✅          |
| [`Zhu`](@ref)                    | Nearest neighbors | `Dataset`  |         ✅          |
| [`ZhuSingh`](@ref)               | Nearest neighbors | `Dataset`  |         ✅          |
| [`Gao`](@ref)                    | Nearest neighbors | `Dataset`  |         ✅          |
| [`Goria`](@ref)                  | Nearest neighbors | `Dataset`  |         ✅          |
| [`Lord`](@ref)                   | Nearest neighbors | `Dataset`  |         ✅          |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors | `Dataset`  |         ✅          |

```@docs
MutualInformationEstimator
KraskovStögbauerGrassberger1
KraskovStögbauerGrassberger2
GaoKannanOhViswanath
GaoOhViswanath
```
