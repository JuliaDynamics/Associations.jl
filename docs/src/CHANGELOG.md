# Changelog for CausalityTools.jl

## Release v0.4.0

### New functionality

- Added [`ExactSimplexIntersectionTest`](@ref) causality test.

## Release v0.3.0

### New functionality

The workflow for estimating causality from time series and dynamical systems has been completely overhauled. New features:

- Common interface to causality testing using the `causality` function and its method. It handles any
    of the following causality tests:
    1. [`CrossMappingTest`](@ref)
    2. [`ConvergentCrossMappingTest`](@ref)
    3. [`JointDistanceDistributionTest`](@ref)
    4. [`JointDistanceDistributionTTest`](@ref)
    5. [`VisitationFrequencyTest`](@ref)
    6. [`TransferOperatorGridTest`](@ref)
    7. [`ApproximateSimplexIntersectionTest`](@ref)

- Integration with [UncertainData.jl](@ref causality_time_series), allowing easy resampling of
    uncertain data to obtain uncertainty estimates on causality statistics. This is done
    by accepting uncertain data as inputs to `causality`.

- Integration with [DynamicalSystems.jl](@ref causality_dynamical_systems). `causality` accepts
    `DiscreteDynamicalSystem`s or `ContinuousDynamicalSystem`s as inputs.

- Syntax overhaul for [low-level estimators](@ref syntax_overview).

- Updated library of [example systems](@ref example_systems).