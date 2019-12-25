# [List of `CausalityTest`s](@id causality_tests_overview)

Most time series causality statistics, although many are so-called non-parameteric,
require the analyst to provide some parameters for their estimation. Outputs 
may be vastly different based on the parameters chosen for a test.

To systematically deal with this fact, we create estimator-specific
[composite types](https://docs.julialang.org/en/v1/manual/types/#Composite-Types-1) for 
the different causality statistics. Every estimator is represented by a compositite types
whose fields are the parameters of the test. A causality test is therefore uniquely defined 
by its parameters. Available causality tests are listed below.

The advantage of this type-based approach is that we can have a unified syntax to 
estimate time series causality statistics, using the [`causality(source::AbstractVector, target::AbstractVector, test::CausalityTest)`](@ref) method.

!!! note
    It is recommended to use the `causality` methods to estimate causality statistics from 
    your time series. However, if you wish to use the low-level methods directly, then these
    are available from the method-specific pages in the meny.

## Distance based tests

- [`CrossMappingTest`](@ref)
- [`ConvergentCrossMappingTest`](@ref)
- [`JointDistanceDistributionTest`](@ref)
- [`JointDistanceDistributionTTest`](@ref)
- [`SMeasureTest`](@ref)

## Entropy based tests

- [`VisitationFrequencyTest`](@ref)
- [`TransferOperatorGridTest`](@ref)
- [`ExactSimplexIntersectionTest`](@ref)
- [`ApproximateSimplexIntersectionTest`](@ref)
- [`NearestNeighbourMITest`](@ref)

## Predictive asymmetry tests

- [`PredictiveAsymmetryTest`](@ref)
- [`NormalisedPredictiveAsymmetryTest`](@ref)

## Applying the tests

All tests can be used to compute causality statistics on empirical 
[scalar valued time series](@ref causality_tests), 
[uncertain time series](@ref causality_uncertain_naiveresampling) and directly 
[from dynamical systems](@ref causality_dynamical_systems) through the `causality` method.

## Schemes for uncertainty handling

If you want to apply a causality test to time series with uncertainties, you can combine any of the
tests above with the provided [uncertainty handling, subsampling and resampling schemes](@ref causality_uncertaindata).
