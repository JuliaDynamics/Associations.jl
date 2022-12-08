# Relative entropy

Relative entropies, or divergences, as they're often called in the literature come in
several forms, depending on what type of entropy they are based on. In CausalityTools.jl,
we indicate the different types of divergences by passing an [`Entropy`](@ref) instance as
the first argument [`entropy_relative`](@ref).

For example, to compute the Shannon relative entropy, do
`entropy_relative(Shannon(), est, x)`, and to compute Rényi relative entropy, do
`entropy_relative(Renyi(; q), est, x)`.

## Implementations

Any of the following estimators can be used with [`entropy_relative`](@ref).

| Estimator                   | Principle         | Input data                       | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) |
| --------------------------- | ----------------- | -------------------------------- | :---------------: | :-------------: | :---------------: |
| [`Wang`](@ref)              | Nearest neighbors | `Dataset`s with equal dimensions |    Continuous     |        x        |         x         |
| [`WangTransformed`](@ref)   | Nearest neighbors | `Dataset`s with equal dimensions |    Continuous     |        x        |         x         |
| [`PoczosSchneiderRE`](@ref) | Nearest neighbors | `Dataset`s with equal dimensions |    Continuous     |   Continuous    |    Continuous     |
|                             |                   |                                  |                   |                 |                   |

## API

```@docs
entropy_relative
RelativeEntropyEstimator
```

## Estimators

```@docs
Wang
WangTransformed
PoczosSchneiderRE
```
