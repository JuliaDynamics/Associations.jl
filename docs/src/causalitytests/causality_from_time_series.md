
# [Causality from time series](@id causality_time_series)

The `causality` function and its methods provide a common interface for testing causal hypotheses.
For analysing time series, all you need to do is provide a `source` and a `target`. Then, choose 
one of the [available causality tests](@ref causality_tests) to quantify the (directional)
dynamical dependence between `source` and `target`.

All high-level causality test are integrated with the uncertainty handling machinery in 
[UncertainData.jl](https://github.com/kahaaga/UncertainData.jl). Any combination of real-valued vectors, `Vector{<:AbstractUncertainValue}`,
or `AbstractUncertainValueDataset` are accepted as inputs to `causality`, making uncertainty
quantification on the causality statistics a breeze.

For more fine grained control over the analysis, check out the [syntax overview](@ref syntax_overview) for low-level estimators.

```@docs
causality
```

## Resampling uncertain time series with constraints

```@docs
causality(source, target,  resampling::ConstrainedResampling, test::CausalityTest)
```

## Resampling schemes

```@docs
ConstrainedResampling
```
