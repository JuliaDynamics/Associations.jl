
# [Entropies API](@id entropies)

The entropies API is defined by

- [`InformationMeasure`](@ref)
- [`entropy`](@ref)
- [`DifferentialInfoEstimator`](@ref)

The entropies API is re-exported from [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl). Why? Continuous/differential versions of many information theoretic
association measures can be written as a function of differential entropy terms, and can
thus be estimated using [`DifferentialInfoEstimator`](@ref)s.

```@docs
ComplexityMeasures.entropy
```

## Definitions

```@docs
InformationMeasure
Shannon
Renyi
Tsallis
Kaniadakis
Curado
StretchedExponential
```

## [`DifferentialInfoEstimator`](@ref)s

CausalityTools.jl reexports [`DifferentialInfoEstimator`](@ref)s from
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).
Why? Any information-based measure that can be written as a function of differential entropies
can be estimated using a [`DifferentialInfoEstimator`](@ref)s. 

```@docs
DifferentialInfoEstimator
```

### Overview

Only estimators compatible with multivariate data are applicable to the multi-argument measures
provided by CausalityTools. Hence, some entropy estimators are missing from the overview
here (see [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl) for
details).

Each [`DifferentialInfoEstimator`](@ref)s uses a specialized technique to approximate relevant
densities/integrals, and is often tailored to one or a few types of generalized entropy.
For example, [`Kraskov`](@ref) estimates the [`Shannon`](@ref) entropy.

| Estimator                    | Principle         | [`Shannon`](@ref) |
| :--------------------------- | :---------------- | :---------------: |
| [`KozachenkoLeonenko`](@ref) | Nearest neighbors |        ✓         |
| [`Kraskov`](@ref)            | Nearest neighbors |        ✓         |
| [`Zhu`](@ref)                | Nearest neighbors |        ✓         |
| [`ZhuSingh`](@ref)           | Nearest neighbors |        ✓         |
| [`Gao`](@ref)                | Nearest neighbors |        ✓         |
| [`Goria`](@ref)              | Nearest neighbors |        ✓         |
| [`Lord`](@ref)               | Nearest neighbors |        ✓         |

### [`Kraskov`](@ref)

```@docs
Kraskov
```

### [`KozachenkoLeonenko`](@ref)

```@docs
KozachenkoLeonenko
```

### [`Zhu`](@ref)

```@docs
Zhu
```

### [`ZhuSingh`](@ref)

```@docs
ZhuSingh
```

### [`Gao`](@ref)

```@docs
Gao
```

### [`Goria`](@ref)

```@docs
Goria
```

### [`Lord`](@ref)

```@docs
Lord
```

## Utilities

```@docs
entropy_maximum
entropy_normalized
```
