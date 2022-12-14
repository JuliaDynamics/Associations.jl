# Divergences (relative entropy)

## API

```@docs
divergence
Divergence
DivergenceDefinition
DivergenceEstimator
```

## Supported divergences

| Estimator                                    | Principle           | Input data                       | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) |
| -------------------------------------------- | ------------------- | -------------------------------- | :---------------: | :-------------: | :---------------: |
| [`RelativeEntropyShannonDifferential`](@ref) | Estimator-dependent | Estimator-dependent              |    Continuous     |        x        |         x         |
| [`RelativeEntropyShannon`](@ref)             | Estimator-dependent | Estimator-dependent              |     Discrete      |        x        |         x         |
| [`RelativeEntropyRenyiDifferential`](@ref)   | Estimator-dependent | Estimator-dependent              |         x         |   Continuous    |         x         |
| [`RelativeEntropyTsallisDifferential`](@ref) | Estimator-dependent | Estimator-dependent              |         x         |        x        |    Continuous     |
| [`RelativeEntropyRenyi`](@ref)               | Estimator-dependent | Estimator-dependent              |         x         |    Discrete     |         x         |
| [`RelativeEntropyTsallis`](@ref)             | Estimator-dependent | Estimator-dependent              |         x         |        x        |     Discrete      |

### Discrete

```@docs
RelativeEntropyShannon
RelativeEntropyRenyi
RelativeEntropyTsallis
TsallisDivergence
RenyiDivergence
ShannonDivergence
```

### Continuous

```@docs
RelativeEntropyShannonDifferential
RelativeEntropyRenyiDifferential
RelativeEntropyTsallisDifferential
TsallisDivergenceDifferential
RenyiDivergenceDifferential
ShannonDivergenceDifferential
```

## Dedicated estimators

Any of the following dedicated estimators can be used with [`divergence`](@ref) in
combination with relevant divergences above.

| Estimator                                    | Principle           | Input data                       | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) |
| -------------------------------------------- | ------------------- | -------------------------------- | :---------------: | :-------------: | :---------------: |
| [`Wang`](@ref)                               | Nearest neighbors   | `Dataset`s with equal dimensions |    Continuous     |        x        |         x         |
| [`WangTransformed`](@ref)                    | Nearest neighbors   | `Dataset`s with equal dimensions |    Continuous     |        x        |         x         |
| [`PoczosSchneiderRE`](@ref)                  | Nearest neighbors   | `Dataset`s with equal dimensions |    Continuous     |   Continuous    |    Continuous     |



```@docs
Wang
WangTransformed
PoczosSchneiderRE
```
