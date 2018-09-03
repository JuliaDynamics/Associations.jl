# CausalityTools

This package provides tools for nonparametric detection of causal relationships between dynamical variables based on observations.

This package brings together the following packages into one environment:


This package provides tools for nonparametric detection of causal relationships between dynamical variables based on observations. The package is equally well-suited for the study of experimental data, and for studying dynamical systems from a more formal context. The workflow integrates nicely  with `DynamicalSystems.jl`.

This package brings together the following packages into one environment:

|Package|Version|Build status|Does what| Provides/builds on | 
|---	|---	|---	|---	|  |
| [`StateSpaceReconstruction.jl`](https://github.com/kahaaga/StateSpaceReconstruction.jl/) | Not registered | [![Build Status](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl.svg?branch=master)](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl) | Flexible state space reconstruction/embeddings. | Instances of some subtype of `AbstractEmbedding` are the foundation for most estimators in the package. |
| [`TimeseriesSurrogates.jl`](https://github.com/kahaaga/TimeseriesSurrogates.jl/) | 0.1.0 | [![Build Status](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl) | Generate surrogate data from time series and orbits. | Null-hypothesis testing. |
| [`TransferEntropy.jl`](https://github.com/kahaaga/TransferEntropy.jl/) | Not registered | [ ![Build Status](https://travis-ci.org/kahaaga/TransferEntropy.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TransferEntropy.jl) | Transfer entropy estimators. | Compute the information flow between time series. All transfer entropy estimators take an embedding as input, discretizes it in some way, then computes the information flow based on that. |
| [`PerronFrobenius.jl`](https://github.com/kahaaga/PerronFrobenius.jl/) | Not registered | [ ![Build Status](https://travis-ci.org/kahaaga/PerronFrobenius.jl.svg?branch=master)](https://travis-ci.org/kahaaga/PerronFrobenius.jl) | Transfer (Perron-Frobenius) operator and invariant measure estimators. | The transfer operator estimators take an embedding as input, discretizes it in some way, then computes the transfer matrix. An invariant measure over the states of the partitioned embedding is computed based on the transfer matrix. |
| [`Simplices.jl`](https://github.com/kahaaga/Simplices.jl/) | Not registered | [![Build Status](https://travis-ci.org/kahaaga/Simplices.jl.svg?branch=master)](https://travis-ci.org/kahaaga/Simplices.jl)|Compute exact simplex intersections in N dimensions. (https://github.com/kahaaga/Simplices.jl/).| Used for the exact triangulation transfer operator estimator in [`Simplices.jl`] |
