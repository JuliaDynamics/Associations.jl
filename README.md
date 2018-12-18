# CausalityTools

| Build status  | Documentation |
| ------------- | ------------- |
| [![Build Status](https://travis-ci.org/kahaaga/CausalityTools.jl.svg?branch=master)](https://travis-ci.org/kahaaga/CausalityTools.jl)  | [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://kahaaga.github.io/CausalityTools.jl/dev)  |

`CausalityTools.jl` provides tools for nonparametric detection of causal relationships between dynamical variables based on time series of observations.

In addition, we provide functions to approximate the transfer operator (Perron-Frobenius operator), and from it, invariant distributions over the discretized state space (embedding).

## What can I do with `CausalityTools`?
The package is equally well-suited both for the study of causal directionality
for experimental data, and for studying dynamical systems from a more formal context. The workflow integrates nicely with `DynamicalSystems.jl`.

Check out the [documentation](https://kahaaga.github.io/CausalityTools.jl/dev) (work in progress!) for more information! Please note that the package is under active development, and that breaking changes may occur until version 1.0 is released.

### Package structure
`CausalityTools.jl` brings together the following packages into one environment:


| package | functionality | version |  build |  
| :---   | :--- |    :---:    |   ---: |  
| [`StateSpaceReconstruction.jl`](https://github.com/kahaaga/StateSpaceReconstruction.jl/) | Fully flexible state space reconstructions (embeddings), partitioning routines (variable-width rectangular, and triangulations), and partition refinement (equal-volume splitting of  simplices). | `0.3.2` | [![Build Status](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl.svg?branch=master)](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl) |
| [`TimeseriesSurrogates.jl`](https://github.com/kahaaga/TimeseriesSurrogates.jl/) | Generate surrogate data from time series. | `0.2.1` | [![Build Status](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl) |
| [`TransferEntropy.jl`](https://github.com/kahaaga/TransferEntropy.jl/) | Transfer entropy estimators. | `0.3.3` | [![Build Status](https://travis-ci.org/kahaaga/TransferEntropy.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TransferEntropy.jl) |  |
| [`PerronFrobenius.jl`](https://github.com/kahaaga/PerronFrobenius.jl/) |  Transfer (Perron-Frobenius) operator estimators. | `0.2.3`  | [![Build Status](https://travis-ci.org/kahaaga/PerronFrobenius.jl.svg?branch=master)](https://travis-ci.org/kahaaga/PerronFrobenius.jl) |
| [`Simplices.jl`](https://github.com/kahaaga/Simplices.jl/) | Exact simplex intersections in N dimensions. | `0.2.2` | [![Build Status](https://travis-ci.org/kahaaga/Simplices.jl.svg?branch=master)](https://travis-ci.org/kahaaga/Simplices.jl) |
| [`CrossMappings.jl`](https://github.com/kahaaga/CrossMappings.jl/) | Exact simplex intersections in N dimensions. | `0.2.3` | [![Build Status](https://travis-ci.org/kahaaga/CrossMappings.jl.svg?branch=master)](https://travis-ci.org/kahaaga/CrossMappings.jl) |


## Wrappers for common use cases
Standard wrappers for the causality detection tools are available for direct application to time series. If you're starting out, these wrappers cover the most common use cases.

For more in-depth analysis, the package comes with state space reconstruction (embedding) and discretization routines, which can also be provided seamlessly to the causality estimators.

*Be careful with using the wrappers: for any real application, you should know what you're doing and utilize the underlying functions which give full control over the embeddings and analysis parameters.*


## Installing
To install the package, run the following lines in the Julia console

```julia
using Pkg
Pkg.add("CausalityTools")
```
