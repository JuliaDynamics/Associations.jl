# CausalityTools.jl

A Julia package for the quantification of invariant measures and the  detecting of information flow or dynamical coupling between time series.

## What can I do with `CausalityTools`?
The package features common causality detection algorithms, such as transfer entropy and related information flows. It is equally well-suited both for the study of experimental data, and for studying dynamical systems from a more formal context. The workflow integrates nicely with `DynamicalSystems.jl`.

Check out the different causality measures in the menu, and feel free to check
out one of the tutorials!


### Package structure
`CausalityTools.jl` build on many separate components, and brings together the following packages into one environment:


| package | version |  build |  
| :---   |    :---:    |   ---: |  
| [`StateSpaceReconstruction.jl`](https://github.com/kahaaga/StateSpaceReconstruction.jl/) | unregistered | [![Build Status](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl.svg?branch=master)](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl) |
| [`TimeseriesSurrogates.jl`](https://github.com/kahaaga/TimeseriesSurrogates.jl/) | `0.1.0` | [![Build Status](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl) |
| [`TransferEntropy.jl`](https://github.com/kahaaga/TransferEntropy.jl/) | unregistered | [![Build Status](https://travis-ci.org/kahaaga/TransferEntropy.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TransferEntropy.jl) |
| [`PerronFrobenius.jl`](https://github.com/kahaaga/PerronFrobenius.jl/) | unregistered  | [![Build Status](https://travis-ci.org/kahaaga/PerronFrobenius.jl.svg?branch=master)](https://travis-ci.org/kahaaga/PerronFrobenius.jl) |
| [`Simplices.jl`](https://github.com/kahaaga/Simplices.jl/) | unregistered  | [![Build Status](https://travis-ci.org/kahaaga/Simplices.jl.svg?branch=master)](https://travis-ci.org/kahaaga/Simplices.jl) |
