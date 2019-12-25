
# [Testing for causality from scalar time series](@id causality_tests)

## [Syntax](@id causality_time_series)

```@docs
causality(source::AbstractVector, target::AbstractVector, test::CausalityTest)
```

Returned values from `causality` depend on the `test` type - see the documentation
for specific tests for details.

### [Uncertainty handling](@id uncertainty_handling)

- All high-level causality tests are integrated with the uncertainty handling 
    machinery in [UncertainData.jl](https://github.com/kahaaga/UncertainData.jl). See the list of 
    [uncertainty handling strategies](@ref causality_uncertaindata) for more details.
- Any combination of real-valued vectors, `Vector{<:AbstractUncertainValue}`, 
    or `AbstractUncertainValueDataset` are accepted as inputs to 
    [`causality`](@ref causality_time_series), making uncertainty quantification on 
    the causality statistics a breeze.

## Why wrapper types?

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
