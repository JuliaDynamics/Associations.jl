# CausalityTools

[![CI](https://github.com/juliadynamics/CausalityTools.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/CausalityTools.jl/actions)
[![](https://img.shields.io/badge/docs-latest_tagged-blue.svg)](https://juliadynamics.github.io/CausalityTools.jl/stable/)
[![](https://img.shields.io/badge/docs-dev_(master)-blue.svg)](https://juliadynamics.github.io/CausalityTools.jl/dev/)
[![codecov](https://codecov.io/gh/JuliaDynamics/CausalityTools.jl/branch/master/graph/badge.svg?token=0b71n6x6AP)](https://codecov.io/gh/JuliaDynamics/CausalityTools.jl)
[![DOI](https://zenodo.org/badge/135443027.svg)](https://zenodo.org/badge/latestdoi/135443027)

CausalityTools.jl is a package for quantifying associations and dynamical coupling
between datasets, independence testing and causal inference.

All further information is provided in the
[documentation](https://juliadynamics.github.io/CausalityTools.jl/dev), which you can either
find online or build locally by running the `docs/make.jl` file.

## Key features

- Association measures from conventional statistics, information theory and dynamical
    systems theory, for example distance correlation, mutual information, transfer entropy,
    convergent cross mapping and a lot more!
- A dedicated API for independence testing, which comes with automatic compatibility with
    every measure-estimator combination you can think of. For example, for conditional
    indepencence testing, we offer the generic [`SurrogateAssociationTest`](@ref)
    (which is fully compatible with [TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl)), and the [`LocalPermutationTest`](@ref).
- A dedicated API for causal network inference based on these measures and independence
    tests.

## Installation

To install the package, run `import Pkg; Pkg.add("CausalityTools")`.
