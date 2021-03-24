# CausalityTools

[![CI](https://github.com/juliadynamics/CausalityTools.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/CausalityTools.jl/actions)
[![](https://img.shields.io/badge/docs-latest_tagged-blue.svg)](https://juliadynamics.github.io/CausalityTools.jl/stable/)
[![](https://img.shields.io/badge/docs-dev_(master)-blue.svg)](https://juliadynamics.github.io/CausalityTools.jl/dev/)
[![codecov](https://codecov.io/gh/JuliaDynamics/CausalityTools.jl/branch/master/graph/badge.svg?token=6XlPGg5nRG)](https://codecov.io/gh/JuliaDynamics/CausalityTools.jl)


`CausalityTools.jl` provides tools for nonparametric detection of causal relationships between dynamical variables based on time series of observations.

Check out the [documentation](https://juliadynamics.github.io/CausalityTools.jl/dev) for more information!

## Key features

Informatic theoretic measures:

- Entropies
- Mutual information
- Transfer entropy
- Predictive asymmetry

Distance-based measures:

- Convergent cross mapping
- S-measure
- Joint distance distribution

<!-->
The package is equally well-suited both for the study of causal directionality
for experimental data, and for studying theoretical systems. Key features include:

- Integration with [UncertainData.jl](https://github.com/kahaaga/UncertainData.jl), which makes
    [handling uncertainties](https://kahaaga.github.io/CausalityTools.jl/dev/causalitytests/causality_from_time_series/) in your causality analyses a breeze.

- Integration with [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl),
    which allows the causality statistics to be 
    [applied to dynamical systems](https://kahaaga.github.io/CausalityTools.jl/dev/causalitytests/causality_from_dynamical_systems/).

- [Library of coupled dynamical systems](https://kahaaga.github.io/CausalityTools.jl/dev/example_systems/example_systems_overview/) 
    for testing algorithm performance.
-->

## Installation

CausalityTools.jl is a registered julia package, you can therefore add the latest tagged release
by running the following lines in the Julia console.

```julia
import Pkg; Pkg.add("CausalityTools")
```

For the latest development version of the package, add the package by referring directly to the GitHub repository.

```julia
import Pkg; Pkg.add("https://github.com/juliadynamics/CausalityTools.jl/")
```

### Fixing SpecialFunction.jl error

During installation, you might get an error related to `SpecialFunctions.jl`. If so, just 
run the following:

```julia
Pkg.build("SpecialFunctions")
```

Then run

```julia
using CausalityTools
```

to verify that everything works.
