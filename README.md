# CausalityTools

This package provides tools for nonparametric detection of causal relationships between dynamical variables based on observations.

## What can I do with `CausalityTools`?
The package implements common causality detection algorithms, such as transfer entropy and related information flows.

It is equally well-suited both for the study of experimental data, and for studying dynamical systems from a more formal context. The workflow integrates nicely with `DynamicalSystems.jl`.

Check out the documentation (coming) for more information!


### Package structure
`CausalityTools.jl` brings together the following packages into one environment:


| package | version |  build |  
| :---   |    :---:    |   ---: |  
| [`StateSpaceReconstruction.jl`](https://github.com/kahaaga/StateSpaceReconstruction.jl/) | unregistered | [![Build Status](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl.svg?branch=master)](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl) |
| [`TimeseriesSurrogates.jl`](https://github.com/kahaaga/TimeseriesSurrogates.jl/) | `0.1.0` | [![Build Status](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl) |
| [`TransferEntropy.jl`](https://github.com/kahaaga/TransferEntropy.jl/) | unregistered | [![Build Status](https://travis-ci.org/kahaaga/TransferEntropy.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TransferEntropy.jl) |
| [`PerronFrobenius.jl`](https://github.com/kahaaga/PerronFrobenius.jl/) | unregistered  | [![Build Status](https://travis-ci.org/kahaaga/PerronFrobenius.jl.svg?branch=master)](https://travis-ci.org/kahaaga/PerronFrobenius.jl) |
| [`Simplices.jl`](https://github.com/kahaaga/Simplices.jl/) | `0.1.0` | [![Build Status](https://travis-ci.org/kahaaga/Simplices.jl.svg?branch=master)](https://travis-ci.org/kahaaga/Simplices.jl) |

Until `PerronFrobenius.jl` and `TransferEntropy.jl` is on METADATA, you can
install these packages by running the following lines in the Julia console:

```julia
Pkg.clone("https://github.com/kahaaga/TransferEntropy.jl")
Pkg.clone("https://github.com/kahaaga/PerronFrobenius.jl")
```
