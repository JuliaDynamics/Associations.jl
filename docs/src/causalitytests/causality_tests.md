
# [Causality tests](@id causality_tests)

The syntax for performing a causality tests is always the same: `causality(source, target, test)`.
Here, `test` is an instance of one of the concrete types listed below. Each of these are containers
for test parameters. Returned values from `causality` depend on the `test` type - see the documentation
for specific tests for details.

## Uncertainty handling

All tests can be used in conjunction with a binning scheme with [`BinnedDataCausalityTest`](@ref) to perform the tests on data with uncertain index (e.g time) values. See also [tutorials](@ref tutorials_causality_from_uncertain_data).

## Distance-based causality tests

```@docs
DistanceBasedCausalityTest
```

## Entropy-based tests

```@docs
EntropyBasedCausalityTest
```

### Transfer entropy tests

```@docs
TransferEntropyCausalityTest
```
