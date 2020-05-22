# [FAQ](@id FAQ)

## I want to see an overview of available causality tests

See the [package overview](@ref package_overview).

## I want to use low-level methods directly

### Transfer entropy

- I want to understand how [marginals are assigned](@ref te_assigning_marginals) during transfer entropy computations.
- [`transferentropy(::Any, ::TEVars, ::RectangularBinning, ::BinningTransferEntropyEstimator)`](@ref) for estimating transfer entropy by discretizing the [state space reconstruction](@ref custom_delay_reconstruction) of your data into a rectangular grid.
- [`transferentropy(::Any, ::TEVars, ::NearestNeighbourMI)`](@ref) for estimating transfer entropy by counting nearest neighbours in [state space reconstruction](@ref custom_delay_reconstruction) of your data.
- [`transferentropy(::Any, ::TEVars, ::RectangularBinning, ::BinningTransferEntropyEstimator)`](@ref) for estimating transfer entropy by discretizing the [state space reconstruction](@ref custom_delay_reconstruction) of your data into a triangular partition.

### Joint distance distribution

- [`joint_distance_distribution(source, target)`](@ref) computes the joint distance distribution over a set of subintervals.
- [`joint_distance_distribution(test::OneSampleTTest, source, target)`](@ref) computes the joint distance distribution over a set of subintervals and returns a t-test for the distribution.

## Cross mappings

- [`crossmap(driver, response)`](@ref) computes cross mapping skills for a given library size.
- [`convergentcrossmap(driver, response, timeseries_lengths)`](@ref) computes cross mapping skills for a given library size, but over multiple time series lengths.

## Regular time series

I have regular time series want to see a tuturial on a specific causality detection method:

- [`SMeasureTest` tutorial](@ref tutorial_SMeasureMest).

## Time series with uncertainties

My data have uncertainties in time indices (and potentially values). I want to use the
built-in uncertainty handling tools to compute a causality statitistic for my data.

### Estimating causality statistics on binned data

- [PredictiveAsymmetryTest](@ref tutorial_PredictiveAsymmetryTest_BinnedMeanResampling).
