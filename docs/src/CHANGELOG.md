# Changelog for CausalityTools.jl

## Release v0.9.2

### Documentation

- Fixed formatting for equations of motion for the Ikeda system.

## Release v0.9.1

### Documentation

- Figures in `PredictiveAsymmetryTest` example are now larger.
- [`NearestNeighbourMITest`](@ref) and [`NormalisedPredictiveAsymmetryTest`](@ref) added to documentation
    for [`EntropyBasedCausalityTest`](@ref).
- Documentation for transfer entropy estimators updated with references to where they are applicable.
- Enlarged some figures.
- Improved tutorial for `SMeasureTest`.

### Non-breaking changes

- Make sure all methods use `source`, `target` and `cond` internally (instead of `x`, `y`, and `z`). 
    These have clearer meanings.

## Release v0.9.0

Release v0.9.0 introduces some new functionality and significant documentation improvements.

### New functionality

#### TransferEntropy.jl

- Added [`NearestNeighbourMI`](@ref) transfer entropy estimator.
- Added [`transferentropy(::Any, ::TEVars, ::NearestNeighbourMI)`](@ref) method.
- Transfer entropy estimators now contain a field `b` which gives the base of the logarithm
    used during transfer entropy computations, and hence dictates the unit of the transfer 
    entropy. By default, `b = 2`, which gives the transfer entropy in bits. The keyword `b` 
    is thus obsolete in all transfer entropy methods that used it before.

#### CausalityTools.jl

- Added [`NearestNeighbourMITest`](@ref) transfer entropy test.
- Added [`NormalisedPredictiveAsymmetryTest`](@ref), for which the predictive asymmetry 
    is normalised to some fraction of the mean of the raw values of the predictive statistic.
- All subtypes of `CausalityTest` are now mutable. This allows adjusting test 
    parameters during sensitivity analyses without creating a new test every time.

#### TimeseriesSurrogates.jl

- Plotting functionality is behind a Requires-block. To use the plotting functionality, you 
    now need to do `using Plots` beforehand.

### Documentation

- Fixed small error in documentation for `DiscreteSystemSetup`.
- Created high-level pages for transfer entropy. Re-did documentation.
- Created high-level pages for cross mappings.

## Release v0.7.1

- Small documentation fix for causality on time series.

## Release v0.7.0

### Breaking changes

- The `ν` parameter for the `ConvergentCrossMappingTest` and `CrossMappingTest` 
    has been changed to `η` to conform with the syntax for transfer entropy tests.

### New features

- All subtypes of `CausalityTest` now have a type parameter `N` indicating
    the number of elements that are returned when applying the tests.

## Release v0.6.1

### Documentation

- Improved documentation on binning schems.
- Improved documentation on delay reconstructions.

## Release v0.6.0

### New features

- Added `SMeasureTest` causality test.
- Added tutorial for `SMeasureTest`.

### Documentation

- Improved documentation.

## Release v0.5.2

### Documentation

- Fixed documentation for `DiscreteSystemSetup`.
- Added section in README on how to deal with SpecialFunction.jl build error.
- Added `Pkg` to dependencies, so that build script executes properly.

## Release v0.5.1

### Tutorials

- Added tutorial for `InterpolateBinTest`.

### Documentation

- Fixed documentation for `DiscreteSystemSetup`.

## Release v0.5.0

- Added `InterpolateBinTest`.
- Added `RandomSequencesTest`.

## Release v0.4.1

### New functionality

- Added `ConstrainedTest`.

### Improvements

- Better documentation for uncertain data causality tests.

## Release v0.4.0

### Bug fixes

- Installing the package using `Pkg.add("CausalityTools")` from the central repository caused an  installation error related to `SpecialFunctions.jl`. Building `SpecialFunctions.jl` is now added to the build step of `CausalityTools`.

### New functionality

- Added [`ExactSimplexIntersectionTest`](@ref) causality test.

## Release v0.3.0

### New functionality

The workflow for estimating causality from time series and dynamical systems has been completely overhauled. New features:

#### Common interface to causality testing using the `causality` function and its methods. 

It handles any of the following causality tests:

1. [`CrossMappingTest`](@ref)

2. [`ConvergentCrossMappingTest`](@ref)

3. [`JointDistanceDistributionTest`](@ref)

4. [`JointDistanceDistributionTTest`](@ref)

5. [`VisitationFrequencyTest`](@ref)

6. [`TransferOperatorGridTest`](@ref)

7. [`ApproximateSimplexIntersectionTest`](@ref)

#### Other features

- Integration with [UncertainData.jl](@ref causality_time_series), allowing easy resampling of
    uncertain data to obtain uncertainty estimates on causality statistics. This is done
    by accepting uncertain data as inputs to `causality`.

- Integration with [DynamicalSystems.jl](@ref causality_dynamical_systems). `causality` accepts
    `DiscreteDynamicalSystem`s or `ContinuousDynamicalSystem`s as inputs.

- Syntax overhaul for [low-level estimators](@ref syntax_overview).

- Updated library of [example systems](@ref example_systems).
