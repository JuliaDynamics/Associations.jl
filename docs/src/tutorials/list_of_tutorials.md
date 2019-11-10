# Tutorial notebooks

The following tutorials are also available as [Juputer notebooks](https://github.com/kahaaga/CausalityTools.jl/tree/master/tutorials).

## Bivariate causality analyses on scalar time series

- [`SMeasureTest` on scalar time series](@ref tutorial_smeasuretest).

## [Causality from uncertain data](@id tutorials_causality_from_uncertain_data)

- [`PredictiveAsymmetryTest` with naive resampling](https://github.com/kahaaga/CausalityTools.jl/blob/master/tutorials/CausalityTools_tutorial_causality_test_on_uncertain_data_predictive_asymmetry.ipynb)

## [Causality from uncertain data, binning](@id tutorials_causality_from_uncertain_data)

Many causality tests require data that are equally spaced in time. If your time series have
timing uncertainties, they can be brought on a common, regularly spaced time grid 
by using a [binned resampling scheme](https://kahaaga.github.io/UncertainData.jl/dev/resampling/resampling_schemes/resampling_with_schemes_uncertain_indexvalue_collections/#binned_resampling_schemes). Then, define an instance of a [`BinnedDataCausalityTest`](@ref) with the desired binning scheme and [causality test](@ref causality_tests) of your choice, and apply the test over one or more draws of the binned dataset. Below are some examples:

- [`PredictiveAsymmetryTest` with `BinnedResampling`](@ref tutorial_PredictiveAsymmetryTest_BinnedResampling) [[notebook](https://github.com/kahaaga/CausalityTools.jl/blob/master/tutorials/CausalityTools_tutorial_BinnedDataCausalityTest_PredictiveAsymmetryTest_BinnedResampling.ipynb)]
- [`PredictiveAsymmetryTest` with `BinnedMeanResampling`](@ref tutorial_PredictiveAsymmetryTest_BinnedMeanResampling) [[notebook](https://github.com/kahaaga/CausalityTools.jl/blob/master/tutorials/CausalityTools_tutorial_BinnedDataCausalityTest_PredictiveAsymmetryTest_BinnedMeanResampling.ipynb)]