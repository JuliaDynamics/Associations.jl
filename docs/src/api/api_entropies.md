
# [Entropies API](@id entropies)

The entropies API is defined by

- [`EntropyDefinition`](@ref)
- [`entropy`](@ref)
- [`DifferentialEntropyEstimator`](@ref)

The entropies API is re-exported from [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl). Why? Continuous/differential versions of many information theoretic
association measures can be written as a function of differential entropy terms, and can
thus be estimated using [`DifferentialEntropyEstimator`](@ref)s.

## Definitions

```@docs
EntropyDefinition
Shannon
Renyi
Tsallis
Kaniadakis
Curado
StretchedExponential
```

### Discrete entropy

```@docs
entropy(::EntropyDefinition, ::ProbabilitiesEstimator, ::Any)
entropy_maximum
entropy_normalized
```

### Differential/continuous entropy

```@docs
DifferentialEntropyEstimator
entropy(::DifferentialEntropyEstimator, ::Any)
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

#### Compatibility table

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
