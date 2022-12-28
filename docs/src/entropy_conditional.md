# Conditional entropy

## Implementations

Any of the following estimators can be used with [`entropy_conditional`](@ref).

| Estimator                                       | Principle           | Input data          | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) |
| ----------------------------------------------- | ------------------- | ------------------- | :---------------: | :-------------: | :---------------: |
| [`ConditionalEntropyShannon`](@ref)             | Estimator-dependent | Estimator-dependent |     Discrete      |        x        |         x         |
| [`ConditionalEntropyRenyi`](@ref)               | Estimator-dependent | Estimator-dependent |         x         |    Discrete     |         x         |
| [`ConditionalEntropyTsallis`](@ref)             | Estimator-dependent | Estimator-dependent |         x         |        x        |     Discrete      |
| [`ConditionalEntropyShannonDifferential`](@ref) | Estimator-dependent | Estimator-dependent |    Continuous     |        x        |         x         |

## API

```@docs
entropy_conditional
ConditionalEntropy
ConditionalEntropyDefinition
ConditionalDifferentialEntropyEstimator
```

## Dedicated estimators

No dedicated estimators of conditional entropy are currently implemented. Use
one of the derived estimators instead.

## Derived estimators

```@docs
ConditionalEntropyShannon
ConditionalEntropyRenyi
ConditionalEntropyTsallis
ConditionalEntropyShannonDifferential
```
