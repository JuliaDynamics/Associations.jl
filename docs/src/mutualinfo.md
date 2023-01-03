# Mutual information

The function that estimates mutual information from data is [`mutualinfo`](@ref).
It can estimate different types of mutual informations from the scientific literature,
and each is represnted here as a subtype of [`MutualInformation`](@ref). Because a
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
MIDefinitionShannonDoubleSum
MIDefinitionShannonH3
```

### Tsallis mutual information

```@docs
MITsallis
MIDefinitionTsallisH3Furuichi
MIDefinitionTsallisH3Martin
```

## [Overview table](@id mutualinfo_overview)

### Continuous/differential mutual information

| Estimator                              | Principle         | [`MIShannon`](@ref) |
| -------------------------------------- | ----------------- | :-----------------: |
| [`KraskovStögbauerGrassberger1`](@ref) | Nearest neighbors |         ✅          |
| [`KraskovStögbauerGrassberger2`](@ref) | Nearest neighbors |         ✅          |
| [`GaoKannanOhViswanath`](@ref)         | Nearest neighbors |         ✅          |
| [`GaoOhViswanath`](@ref)               | Nearest neighbors |         ✅          |

Continuous/differential mutual information may also be estimated using any of our
differential entropy estimators that support multivariate input data.

| Estimator                        | Principle         | Input data | [`MIShannon`](@ref) |
| -------------------------------- | ----------------- | ---------- | :-----------------: |
| [`Kraskov`](@ref)                | Nearest neighbors | `Dataset`  |         ✅          |
| [`Zhu`](@ref)                    | Nearest neighbors | `Dataset`  |         ✅          |
| [`ZhuSingh`](@ref)               | Nearest neighbors | `Dataset`  |         ✅          |
| [`Gao`](@ref)                    | Nearest neighbors | `Dataset`  |         ✅          |
| [`Goria`](@ref)                  | Nearest neighbors | `Dataset`  |         ✅          |
| [`Lord`](@ref)                   | Nearest neighbors | `Dataset`  |         ✅          |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors | `Dataset`  |         ✅          |

### Discrete estimators

Discrete mutual information may also be estimated using any of the following
probability estimators (the table scrolls sideways). Instead 

For some of these methods to work, specialized treatment of the input data is required.
The idea is to 

| Estimator                    | Principle                   |     Input data     | [`MIDefinitionShannonH3`](@ref) | [`MIDefinitionTsallisH3Furuichi`](@ref) | [`MIDefinitionTsallisH3Martin`](@ref) |
| ---------------------------- | --------------------------- | :----------------: | :-----------------------------: | :-------------------------------------: | :-----------------------------------: |
| [`CountOccurrences`](@ref)   | Frequencies                 | `Vector`/`Dataset` |               ✅                |                   ✅                    |                  ✅                   |
| [`ValueHistogram`](@ref)     | Binning (histogram)         | `Vector`/`Dataset` |               ✅                |                   ✅                    |                  ✅                   |
| [`TransferOperator`](@ref)   | Binning (transfer operator) | `Vector`/`Dataset` |               ✅                |                   ✅                    |                  ✅                   |
| [`SymbolicPermuation`](@ref) | Ordinal patterns            |      `Vector`      |               ✅                |                   ✅                    |                  ✅                   |
| [`Dispersion`](@ref)         | Dispersion patterns         |      `Vector`      |               ✅                |                   ✅                    |                  ✅                   |

## Estimators

```@docs
MutualInformationEstimator
KraskovStögbauerGrassberger1
KraskovStögbauerGrassberger2
GaoKannanOhViswanath
GaoOhViswanath
```
