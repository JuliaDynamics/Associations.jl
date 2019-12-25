
# [Testing for causality from scalar time series](@id causality_tests)

## [Syntax](@id causality_time_series)

```@docs
causality(source::AbstractVector, target::AbstractVector, test::CausalityTest)
```

Returned values from `causality` depend on the `test` type - see the documentation
for [specific tests](@ref causality_tests_overview) for details.

## [Uncertainty handling](@id uncertainty_handling)

- All high-level causality tests are integrated with the uncertainty handling 
    machinery in [UncertainData.jl](https://github.com/kahaaga/UncertainData.jl). See the list of 
    [uncertainty handling strategies](@ref causality_uncertaindata) for more details.
- Any combination of real-valued vectors, `Vector{<:AbstractUncertainValue}`, 
    or `AbstractUncertainValueDataset` are accepted as inputs to 
    [`causality`](@ref causality_time_series), making uncertainty quantification on 
    the causality statistics a breeze.
