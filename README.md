# CausalityTools

[![CI](https://github.com/juliadynamics/TransferEntropy.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/TransferEntropy.jl/actions)
[![](https://img.shields.io/badge/docs-latest_tagged-blue.svg)](https://juliadynamics.github.io/TransferEntropy.jl/stable/)
[![](https://img.shields.io/badge/docs-dev_(master)-blue.svg)](https://juliadynamics.github.io/TransferEntropy.jl/dev/)
[![codecov](https://codecov.io/gh/JuliaDynamics/TransferEntropy.jl/branch/master/graph/badge.svg?token=6XlPGg5nRG)](https://codecov.io/gh/JuliaDynamics/TransferEntropy.jl)

`CausalityTools.jl` provides methods for causal inference and detection of dynamical coupling based on time series.

Check out the [documentation](https://juliadynamics.github.io/CausalityTools.jl/dev) for more information!

## Key tools

- A easy-to-use framework for estimating information theoretic measures, such as transfer entropy, predictive asymmetry.
- Generalized entropy and mutual information.
- Convergent cross mapping, S-measure and joint distance distribution.
- Surrogate data generation.

## Installation

CausalityTools.jl is a registered julia package, you can therefore add the latest tagged release
by running the following lines in the Julia console.

```julia
import Pkg; Pkg.add("CausalityTools")
```

For the latest development version of the package, add the package by referring directly to the GitHub repository.

```julia
import Pkg; Pkg.add(url="https://github.com/juliadynamics/CausalityTools.jl/", rev="master")
```
