```@meta
DocTestSetup = quote
    using CausalityTools
    using Plots
end
```

```@setup overall_setup
using CausalityTools
using Plots
```

# CausalityTools.jl

`CausalityTools` is a Julia package providing algorithms for detecting causal relations in complex systems based on time series data.

## Goals

- **Provide a comprehensive, [easy-to-use framework](@ref causality_time_series)** for the detection of directional causal influences in complex dynamical systems from time series.

- **Functional and efficient [implementations](@ref syntax_overview)** of causality detection algorithms, with thorough documentation and references to primary literature.

- [**Integration with UncertainDatasets.jl**](@ref causality_time_series), which greatly simplifies working with uncertain data.

- [**Integration with DynamicalSystems.jl**](@ref causality_dynamical_systems), for quick analysis of time series from systems where the governing equations are known.

- **[Library of example dynamical systems](@ref example_systems)** for testing algorithm performance.

- **[Surrogate data methods](@ref surrogate_methods)** for null-hypothesis testing. In the future,
    surrogate methods will be provided as part of resampling schemes.

- **Worked examples** for the algorithms.

## Status

The package and documentation is under active development. **Breaking changes may occur in CausalityTools and its dependencies until the 1.0 release**.

### Package structure

`CausalityTools.jl` brings together the following packages into one environment:

| package | functionality | version |  build |  
| :---   | :--- |    :---:    |   ---: |  
| [`CausalityToolsBase.jl`](https://github.com/kahaaga/CausalityToolsBase.jl/) | Basic functionality for the CausalityTools ecosystem. | `0.6.0` | [![Build Status](https://travis-ci.com/kahaaga/CausalityToolsBase.jl.svg?branch=master)](https://travis-ci.com/kahaaga/CausalityToolsBase.jl) |
| [`StateSpaceReconstruction.jl`](https://github.com/kahaaga/StateSpaceReconstruction.jl/) | Fully flexible state space reconstructions (embeddings), partitioning routines (variable-width rectangular, and triangulations), and partition refinement (equal-volume splitting of  simplices). | `0.4.2` | [![Build Status](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl.svg?branch=master)](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl) |
| [`TimeseriesSurrogates.jl`](https://github.com/kahaaga/TimeseriesSurrogates.jl/) | Generate surrogate data from time series. | `0.3.1` | [![Build Status](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TimeseriesSurrogates.jl) |
| [`TransferEntropy.jl`](https://github.com/kahaaga/TransferEntropy.jl/) | Transfer entropy estimators. | `0.4.3` | [![Build Status](https://travis-ci.org/kahaaga/TransferEntropy.jl.svg?branch=master)](https://travis-ci.org/kahaaga/TransferEntropy.jl) |  |
| [`PerronFrobenius.jl`](https://github.com/kahaaga/PerronFrobenius.jl/) |  Transfer (Perron-Frobenius) operator estimators. | `0.6.0`  | [![Build Status](https://travis-ci.org/kahaaga/PerronFrobenius.jl.svg?branch=master)](https://travis-ci.org/kahaaga/PerronFrobenius.jl) |
| [`Simplices.jl`](https://github.com/kahaaga/Simplices.jl/) | Exact simplex intersections in N dimensions. | `0.4.1` | [![Build Status](https://travis-ci.org/kahaaga/Simplices.jl.svg?branch=master)](https://travis-ci.org/kahaaga/Simplices.jl) |
| [`CrossMappings.jl`](https://github.com/kahaaga/CrossMappings.jl/) | Exact simplex intersections in N dimensions. | `0.3.4` | [![Build Status](https://travis-ci.org/kahaaga/CrossMappings.jl.svg?branch=master)](https://travis-ci.org/kahaaga/CrossMappings.jl) |

## Contributors

- Kristian Agasøster Haaga ([@kahaaga](https://github.com/kahaaga))
- David Diego ([@susydiegolas](https://github.com/susydiegolas))
- Tor Einar Møller ([@tormolle](https://github.com/tormolle))

## Related software

- [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl) provides a range of tools for exploring nonlinear dynamics and chaos, both for synthetic and observed systems. We provide seamless interaction with `Dataset` outputs from DynamicalSystems.  Most of our example systems are also implemented as `DiscreteDynamicalSystem`s or `ContinuousDynamicalSystems` from DynamicalSystems.
