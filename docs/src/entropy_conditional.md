# Conditional entropy

Conditional entropies appear in many forms in the literature, depending on what type of 
entropy they are based on. In CausalityTools.jl, we indicate the different types of divergences by
passing an [`Entropy`](@ref) instance as the first argument [`entropy_relative`](@ref).

For example, to compute the Shannon conditional entropy, do
`entropy_conditional(Shannon(), est, x)`, and to compute Rényi conditional entropy, do 
`entropy_conditional(Renyi(; q), est, x)`.

## Implementations

Any of the following estimators can be used with [`entropy_cross`](@ref).

| Estimator            | Principle | Input data          | [`Shannon`](@ref) |   [`Renyi`](@ref)    | [`Tsallis`](@ref) |
| -------------------- | --------- | ------------------- | :---------------: | :------------------: | :---------------: |
| [`FuruichiCD`](@ref) |           | Estimator-dependent |         x         |          x           |     Discrete      |
| [`ShannonCD`](@ref)  |           | Estimator-dependent |         x         | Discrete, continuous |         x         |
| [`Jizba`](@ref)      |           | Estimator-dependent |         x         |       Discrete       |         x         |

## API

```@docs
entropy_cross
ConditionalEntropyEstimator
```

## Estimators

```@docs
FuruichiCD
Jizba
ShannonCD
```
