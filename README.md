# Associations

[![CI](https://github.com/juliadynamics/Associations.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/Associations.jl/actions)
[![](https://img.shields.io/badge/docs-latest_tagged-blue.svg)](https://juliadynamics.github.io/Associations.jl/stable/)
[![](https://img.shields.io/badge/docs-dev_(main)-blue.svg)](https://juliadynamics.github.io/Associations.jl/dev/)
[[![codecov](https://codecov.io/gh/JuliaDynamics/Associations.jl/branch/master/graph/badge.svg?token=0b71n6x6AP)](https://codecov.io/gh/JuliaDynamics/Associations.jl)](https://img.shields.io/codecov/c/github/juliadynamics/Associations.jl/main)
[![DOI](https://zenodo.org/badge/135443027.svg)](https://zenodo.org/badge/latestdoi/135443027)

Associations.jl is a package for quantifying associations, independence testing and causal inference.

All further information is provided in the
[documentation](https://juliadynamics.github.io/Associations.jl/dev), which you can either
find online or build locally by running the `docs/make.jl` file.

## Key features

- **Association API**: includes measures and their estimators for pairwise, conditional and other forms of 
    association from conventional statistics, from dynamical systems theory, and from information theory: partial correlation, distance correlation, (conditional) mutual information, transfer entropy, convergent cross mapping and a lot more!
- **Independence testing API**, which is automatically compatible with
    every association measure estimator implemented in the package. 
- **Causal (network) inference API** integrating the association measures and independence testing framework.

## Addititional features

Extending on features from [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl),
we also offer 

- Discretization API for multiple (multivariate) input datasets.
- Multivariate counting and probability estimation API.
- Multivariate information measure API

## Installation

To install the package, run `import Pkg; Pkg.add("Associations")`.

*Previously, this package was called CausalityTools.jl*.