# Conditional mutual information

Conditional mutual informations come in several forms in the literature, depending on what
type of entropy they are based on. In CausalityTools.jl, we indicate the different types
of conditional mutual informations by passing an [`Entropy`](@ref) instance as the first
argument to [`mutualinfo`](@ref).

For example, to compute the Shannon conditional mutual information, do
`condmutualinfo(Shannon(), est, x, y)`, to compute Tsallis conditional mutual information, do
`condmutualinfo(Tsallis(; q), est, x, y)`.

## Implementations

Any of the following estimators can be used with [`condmutualinfo`](@ref).

| Estimator                           | Principle         | Input data         |  [`Shannon`](@ref)  | [`Renyi`](@ref) | [`Tsallis`](@ref) |
| ----------------------------------- | ----------------- | ------------------ | :-----------------: | :-------------: | :---------------: |
| [`FrenzelPompeVelmejkaPalus`](@ref) | Nearest neighbors | `Vector`,`Dataset` |     Continuous      |        x        |         x         |
| [`Rahimzamani`](@ref)               | Nearest neighbors | `Vector`,`Dataset` | Continuous/discrete |        s        |         x         |
| [`PoczosSchneiderCMI`](@ref)        | Nearest neighbors | `Vector`,`Dataset` |          x          |   Continuous    |         x         |

## API

```@docs
condmutualinfo
ConditionalMutualInformationEstimator
```

## Dedicated estimators

```@docs
Rahimzamani
FrenzelPompeVelmejkaPalus
PoczosSchneiderCMI
```
