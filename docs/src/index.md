# CausalityTools.jl

`CausalityTools` is a Julia package that provides algorithms for *detecting
dynamical influences* and *causal inference* based on time series data, and other
commonly used measures of dependence.

!!! info
    You are reading the development version of the documentation of
    CausalityTools.jl that will become version 2.0.

## Information measures

Information measures are build on probabilities and/or entropies.
CausalityTools.jl has a modular design, meaning that any "high-level" information measures may be expressed in terms of either probabilities or entropies.

- Information measures are computed in their discrete form by using
    [`ProbabilitiesEstimator`](@ref).
- Information measures are computed in their differential/continuous
    form by using dedicated estimators, such as
    [`DifferentialEntropyEstimator`](@ref), or [`MutualInformationEstimator`](@ref).

## Input data

```@docs
Dataset
```
