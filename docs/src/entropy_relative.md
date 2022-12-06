# Relative entropy

Relative entropies, or divergences, as they're often called in the literature come in several forms. In CausalityTools.jl, we indicate the different types of divergences by passing
an [`Entropy`](@ref) instance as the first argument [`entropy_relative`](@ref). For example, to compute the Shannon relative entropy, pass `Shannon()` as the first argument, and to compute Rényi relative entropy, pass `Renyi(; q)` as the first argument.

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

## Implementations

| Estimator                   | Principle         | Input data                       | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) |
| --------------------------- | ----------------- | -------------------------------- | :---------------: | :-------------: | :---------------: |
| [`Wang`](@ref)              | Nearest neighbors | `Dataset`s with equal dimensions |     Continuous      |        x        |         x         |
| [`WangTransformed`](@ref)   | Nearest neighbors | `Dataset`s with equal dimensions |     Continuous      |        x        |         x         |
| [`PoczosSchneiderRE`](@ref) | Nearest neighbors | `Dataset`s with equal dimensions |    Continuous     |   Continuous    |    Continuous     |
|                             |                   |                                  |                   |                 |                   |
