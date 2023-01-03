# [Entropies](@id entropies)

## Entropies API

The entropies API is defined by

- [`EntropyDefinition`](@ref)
- [`entropy`](@ref)
- [`DifferentialEntropyEstimator`](@ref)

Please be sure you have read the [Terminology](@ref) section before going through the API here, to have a good idea of the different "flavors" of entropies and how they all come together over the common interface of the [`entropy`](@ref) function.

## Entropy definitions

```@docs
EntropyDefinition
Shannon
Renyi
Tsallis
Kaniadakis
Curado
StretchedExponential
```

## Discrete entropy

```@docs
entropy(::EntropyDefinition, ::ProbabilitiesEstimator, ::Any)
entropy_maximum
entropy_normalized
```

## Differential entropy

```@docs
entropy(::EntropyDefinition, ::DifferentialEntropyEstimator, ::Any)
```

### Table of differential entropy estimators

The following estimators are *differential* entropy estimators, and can also be used
with [`entropy`](@ref).

Each [`DifferentialEntropyEstimator`](@ref)s uses a specialized technique to approximate relevant
densities/integrals, and is often tailored to one or a few types of generalized entropy.
For example, [`Kraskov`](@ref) estimates the [`Shannon`](@ref) entropy.

| Estimator                    | Principle         | Input data | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) | [`Kaniadakis`](@ref) | [`Curado`](@ref) | [`StretchedExponential`](@ref) |
| :--------------------------- | :---------------- | :--------- | :---------------: | :-------------: | :---------------: | :------------------: | :--------------: | :----------------------------: |
| [`KozachenkoLeonenko`](@ref) | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Kraskov`](@ref)            | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Zhu`](@ref)                | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`ZhuSingh`](@ref)           | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Gao`](@ref)                | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Goria`](@ref)              | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Lord`](@ref)               | Nearest neighbors | `Dataset`  |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Vasicek`](@ref)            | Order statistics  | `Vector`   |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Ebrahimi`](@ref)           | Order statistics  | `Vector`   |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`Correa`](@ref)             | Order statistics  | `Vector`   |        ✓         |        x        |         x         |          x           |        x         |               x                |
| [`AlizadehArghami`](@ref)    | Order statistics  | `Vector`   |        ✓         |        x        |         x         |          x           |        x         |               x                |

```@docs
DifferentialEntropyEstimator
```

```@docs
Kraskov
KozachenkoLeonenko
Zhu
ZhuSingh
Gao
Goria
Lord
Vasicek
AlizadehArghami
Ebrahimi
Correa
```
