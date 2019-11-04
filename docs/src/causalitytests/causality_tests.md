
# [Causality tests](@id causality_tests)

`causality` accepts many different causality tests. The returned values will vary depending 
on the type of test, but the syntax for performing a test is always the same.

All tests can be used in conjunction with a binning scheme with [`BinnedDataCausalityTest`](@ref) to perform the tests on data with uncertain index (e.g time) values. See also [tutorials](@ref tutorials_causality_from_uncertain_data).

## Distance based causality tests

- [`CrossMappingTest`](@ref)

- [`ConvergentCrossMappingTest`](@ref)

- [`JointDistanceDistributionTest`](@ref)

- [`JointDistanceDistributionTTest`](@ref)

## Entropy based causality tests

- [`VisitationFrequencyTest`](@ref)

- [`TransferOperatorGridTest`](@ref)

- [`ApproximateSimplexIntersectionTest`](@ref)

- [`ExactSimplexIntersectionTest`](@ref)


## Predictive asymmetry tests

- [`PredictiveAsymmetryTest`](@ref)
