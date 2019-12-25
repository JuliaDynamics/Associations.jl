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

## Introduction

`CausalityTools` is a Julia package providing algorithms for detecting causal relations in complex systems based on time series data.

## Getting started

Causality statistics and their different methods/estimators are described on method-specific pages. These can be accessed from the 'Package ecosystem' tab in the top menu.

- [Predictive asymmetry](@ref predictive_asymmetry_overview)
- [Transfer entropy](@ref transferentropy_overview)
- [Cross mapping](@ref crossmapping_overview)
- [S-Measure](@ref SMeasureTest_overview)
- [Joint distance distribution](@ref joint_distance_distribution_overview)

For more flexibility, however, check out the [`causality` interface](@ref causality_tests_overview),
which provides a unified syntax for applying all the tests, even on uncertain data!

## Contents

The `CausalityTools` module provides set of tools for time series causality testing. Much of the machinery required for these to work are organised in different packages/submodules. Simpler algorithms are implemented directly as submodules in `CausalityTools`.  Relevant types and methods are all re-exported by `CausalityTools`. The functionality and purpose of these modules are summarised below.

### [CausalityTests](@ref causality_tests_overview)

Provides parameter wrapper types for the different causality statistics. This gives unified syntax to evaluate different statistics on empirical [scalar valued time series](@ref causality_tests), 
[uncertain time series](@ref causality_uncertain_naiveresampling) and directly 
[from dynamical systems](@ref causality_dynamical_systems) through the `causality` method.

For time series, the syntax is [`causality(source, target, test::CausalityTest)`](@ref causality_time_series). Here, `source` and `target` are your time series (may also be uncertain data, see below) and `test` is an instance of any of the following causality test parameter types.

#### Transfer entropy tests

- [`VisitationFrequencyTest`](@ref)
- [`TransferOperatorGridTest`](@ref)
- [`ExactSimplexIntersectionTest`](@ref)
- [`ApproximateSimplexIntersectionTest`](@ref)
- [`NearestNeighbourMITest`](@ref)

#### Predictive asymmetry tests

- [`PredictiveAsymmetryTest`](@ref)
- [`NormalisedPredictiveAsymmetryTest`](@ref)

#### Distance based tests

- [`CrossMappingTest`](@ref)
- [`ConvergentCrossMappingTest`](@ref)
- [`JointDistanceDistributionTest`](@ref)
- [`JointDistanceDistributionTTest`](@ref)
- [`SMeasureTest`](@ref)

For those interested, low-level methods are available individual method 
pages under `Package ecosystem` in the menu. These also provide a bit more 
information on each of the causality statistics.

### [CausalityToolsBase](@ref custom_delay_reconstructions)

- [Generalised delay reconstructions](@ref custom_delay_reconstruction) designed to assist with transfer entropy computations and related methods.
- [Discretization routines](@ref discretization) for partitioning state space reconstructions, or multidimensional datasets in general.

### [Simplices](@ref exact_simplex_intersection)

- Computing exact intersections between n-dimensional simplices. Used internally in the [transfer operators](@ref transferoperator_estimation) estimators.

### [PerronFrobenius](@ref transferoperator_estimation)

- [Transfer operator approximation](@ref transferoperator_estimation) over [rectangular](@ref transferoperator_estimation_rectangular) or [triangulated](@ref transferoperator_estimation_triang) partitions.
- [Estimating invariant measures](@ref invariant_measure_estimation) over [rectangular](@ref invariant_measure_estimation_rectangular) or [triangulated](@ref invariant_measure_estimation_triang) partitions, by using the transfer operator approximations.

### [CrossMappings](@ref ccm)

- Implementations of the (convergent) cross mapping routines from Sugihara et al. (2012).

### [Joint distance distribution](@ref joint_distance_distribution_overview)

- Implementation of the joint distance distribution algorithm for detecting directional causality.

### [SMeasure](@ref Smeasure_overview)

- Implementation of the S-measure algorithm for detecting directional causality.

### [TransferEntropy](@ref te_estimation_procedure)

- Defines transfer entropy [estimator types](@ref te_estimators).
- Convenient [delay reconstructions](@ref te_embed_funcs) specific to transfer entropy estimation.

### [TimeseriesSurrogates](@ref surrogate_methods)

`TimeseriesSurrogates` provides a set of methods to generate surrogate time series. It provides 
functionality useful outside causality testing, so is organised in an independent package which 
also maintains [its own documentation](https://kahaaga.github.io/TimeseriesSurrogates.jl/latest/).

Selected features:

- Random surrogates
- Fourier transform amplitude surrogates
- Fourier transform phase surrogates
- Amplitude-adjusted Fourier transform surrogates.
- Visualisation tools making it easy to visualise surrogate data for your time series.

### [UncertainData](@ref uncertainty_handling)

[`UncertainData`](https://github.com/kahaaga/UncertainData.jl) is a package for generic 
uncertain data handling. It was originally developed for handling uncertainties in proxy 
data when doing causal analyses, but provides functionality useful elsewhere, so is kept 
in a separate package.

- Defines `AbstractUncertainValue` and its concrete subtypes, which represent 
    uncertainties of different types, such as theoretical distributions (e.g. normal 
    distributions with mean and standard deviation), kernel density estimated 
    distributions, fitted distributions, and weighted populations (population members 
    may themself be uncertain values.

- Provides the `UncertainIndexValueDataset` type which allows you to define time series 
    whose data points may have uncertainties both in time indices and values. Different
    uncertainty types can be mixed seamlessly.
- Has routines for resampling, binning and interpolating `UncertainIndexValueDataset`s.

`CausalityTools` is integrated with [`UncertainData`](https://github.com/kahaaga/UncertainData.jl) (see below) in the following way:

- Any combination of  real-valued vectors, `Vector{<:AbstractUncertainValue}`, 
    or `AbstractUncertainValueDataset` are accepted as inputs to 
    [`causality`](@ref causality_time_series), making uncertainty quantification on 
    the causality statistics a breeze.
- Causality test parameter types can be [combined](@ref causality_uncertaindata) 
    with resampling routines from [UncertainData.jl](https://github.com/kahaaga/UncertainData.jl) to 
    provide meta-tests which employs more advanced [resampling schemes](@ref causality_uncertaindata) to compute the causality statistics from your time series.

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

It is almost certainly guaranteed that the high level [`causality`](@ref causality_time_series) methods and [`CausalityTest`] types will persist in their current form. Low-level methods might change slightly. Using the [`causality`](@ref causality_time_series) methods to compute various statistics for your time series is therefore more future-proof than directly applying the low-level methods.

## Installation

`CausalityTools` is a registrered Julia package. You can install it by running the following in a Julia console.

```julia
using Pkg
Pkg.add("CausalityTools")
```

You can also enter package mode from the console by typing `]`, then `add CausalityTools`.

## Main contributors

- Kristian Agasøster Haaga ([@kahaaga](https://github.com/kahaaga))
- David Diego ([@susydiegolas](https://github.com/susydiegolas))
- Tor Einar Møller ([@tormolle](https://github.com/tormolle))

## Related software

- [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl) provides a range of tools for exploring nonlinear dynamics and chaos, both for synthetic and observed systems. We provide seamless interaction with `Dataset` outputs from DynamicalSystems.  Most of our example systems are also implemented as `DiscreteDynamicalSystem`s or `ContinuousDynamicalSystems` from DynamicalSystems.
