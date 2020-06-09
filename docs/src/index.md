# Overview

CausalityTools.jl is a package that provides algorithms for causal inference from time series.

CausalityTools.jl is part of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/), a GitHub organization dedicated to creating high quality scientific software.

The package is maintained by the [Earth System Evolution](https://www.earthsystemevolution.com/) group at the 
University of Bergen, Norway.

## Installation 

CausalityTools.jl is a registered Julia package. Install the latest stable version by running the following in the Julia repl:

```julia
using Pkg; Pkg.add("CausalityTools")
```

For the latest development version, you can also add the master branch directly (but on this branch, breaking changes may occur).

```julia 
using Pkg; Pkg.add("https://github.com/JuliaDynamics/CausalityTools.jl")
```


## Utilities 

Utility methods are also provided for:

- transfer operator approximation 
- invariant distribution estimation
- time series surrogate generation.

## Citing

If using any of the methods for scientific research, please cite the original research referenced in the docstrings. 
Some implementations may be based on more than one paper.