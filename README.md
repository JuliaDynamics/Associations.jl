# CausalityTools

| Build status  | Documentation |
| ------------- | ------------- |
| [![Build Status](https://travis-ci.org/kahaaga/CausalityTools.jl.svg?branch=master)](https://travis-ci.org/kahaaga/CausalityTools.jl)  | [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://kahaaga.github.io/CausalityTools.jl/dev)  |

`CausalityTools.jl` provides tools for nonparametric detection of causal relationships between dynamical variables based on time series of observations.

Check out the [documentation](https://kahaaga.github.io/CausalityTools.jl/dev) for more information!

## Key features

The package is equally well-suited both for the study of causal directionality
for experimental data, and for studying theoretical systems. Key features include:

- integration with [UncertainData.jl](https://github.com/kahaaga/UncertainData.jl), which makes
    [handling uncertainties](https://kahaaga.github.io/CausalityTools.jl/dev/causalitytests/causality_from_time_series/) in your causality analyses a breeze.

- integration with [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl),
    which allows the causality statitics to be 
    [applied to dynamical systems](https://kahaaga.github.io/CausalityTools.jl/dev/causalitytests/causality_from_dynamical_systems/).

- [Library of coupled dynamical systems](https://kahaaga.github.io/CausalityTools.jl/dev/example_systems/example_systems_overview/) 
    for testing algorithm performance.

## Installation

CausalityTools.jl is a registered julia package, you can therefore add the latest tagged release
by running the following lines in the Julia console.

```julia
import Pkg; Pkg.add("CausalityTools")
```

For the latest development version of the package, add the package by referring directly to the GitHub repository.

```julia
import Pkg; Pkg.add("https://github.com/kahaaga/CausalityTools.jl/")
```
