# [List of `CausalityTest`s](@id causality_tests_overview)

Most time series causality statistics, although many are so-called non-parameteric,
require the analyst to provide some parameters for their estimation. Outputs 
may be vastly different based on the parameters chosen for a test.

To systematically deal with this fact, we create estimator-specific
[composite types](https://docs.julialang.org/en/v1/manual/types/#Composite-Types-1) for 
the different causality statistics. Every estimator is represented by a composite types
whose fields are the parameters of the test. A causality test is therefore uniquely defined 
by its parameters. Available causality tests are listed below.

The advantage of this type-based approach is that we can have a unified syntax to 
estimate time series causality statistics, using the [`causality(source::AbstractVector, target::AbstractVector, test::CausalityTest)`](@ref) method, which is all a beginner-user
needs to learn. The output of a call to `causality` depends on which `test` is provided, 
but the input format is always the same: your time series, plus a causality test.

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

Low-level causality statistics and their different methods/estimators are described on the following pages:

- [Predictive asymmetry](@ref predictive_asymmetry_overview)
- [Transfer entropy](@ref transferentropy_overview)
- [Cross mapping](@ref crossmapping_overview)
- [S-Measure](@ref Smeasure_overview)
- [Joint distance distribution](@ref joint_distance_distribution_overview)

## Schemes for uncertainty handling

If you want to apply a causality test to time series with uncertainties, you can combine any of the
tests above with the provided [uncertainty handling, subsampling and resampling schemes](@ref causality_uncertaindata).

## Why wrapper types for the parameters?

There may be many ways of estimating a particular causality statistic. For example, 
there are many  transfer entropy estimators, for example 

- [`VisitationFrequency`](@ref), and
- [`NearestNeigboursMI`](@ref)

Their corresponding  `transferentropy` methods accept different inputs:

- [`transferentropy(pts, vars::TEVars, ϵ::RectangularBinning, estimator::BinningTransferEntropyEstimator)`](@ref)
- [`transferentropy(pts, vars::TEVars, estimator::NearestNeighbourMI)`](@ref)

Wouldn't it be nice to be able to use the same syntax for both estimators? Or for five more estimators 
that quantify causality in different manners? Instead of 
providing test parameters directly to the methods (e.g. [`transferentropy`](@ref)) directly, 
we can create a test instance that contains all the test parameters.

- `test_vf = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = 1:5)`
- `test_nn = NearestNeighbourMITest(ηs = 1:5)`

The `causality` function and its methods provide a common interface for testing causal hypotheses.
For analysing time series, all you need to do is provide a `source` and a `target`. Then, choose 
one of the [available causality tests](@ref causality_tests) to quantify the (directional)
dynamical dependence between `source` and `target`. If your data are uncertain, you can 
use [`resampling schemes`](@ref causality_uncertaindata) to deal with that.
