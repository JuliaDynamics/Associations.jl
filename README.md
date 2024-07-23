# CausalityTools

[![CI](https://github.com/juliadynamics/CausalityTools.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/CausalityTools.jl/actions)
[![](https://img.shields.io/badge/docs-latest_tagged-blue.svg)](https://juliadynamics.github.io/CausalityTools.jl/stable/)
[![](https://img.shields.io/badge/docs-dev_(master)-blue.svg)](https://juliadynamics.github.io/CausalityTools.jl/dev/)
[![codecov](https://codecov.io/gh/JuliaDynamics/CausalityTools.jl/branch/master/graph/badge.svg?token=0b71n6x6AP)](https://codecov.io/gh/JuliaDynamics/CausalityTools.jl)
[![DOI](https://zenodo.org/badge/135443027.svg)](https://zenodo.org/badge/latestdoi/135443027)

CausalityTools.jl is a package for quantifying associations, independence testing and causal inference.

All further information is provided in the
[documentation](https://juliadynamics.github.io/CausalityTools.jl/dev), which you can either
find online or build locally by running the `docs/make.jl` file.

## Key features

- **Association API** including measures and their estimators for pairwise, conditional and other forms of association from conventional statistics, from dynamical
    systems theory, and from information theory: partial correlation, distance correlation, (conditional) mutual information, transfer entropy, convergent cross mapping and a lot more!
- **Independence testing API**, which is automatically compatible with
    every association measure estimator implemented in the package. 
- **Causal (network) inference API** integrating the association measures and independence testing framework.

## Addititional features

- Multivariate probabilities estimation, extending the API from 
    [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).


## Installation

To install the package, run `import Pkg; Pkg.add("CausalityTools")`.
