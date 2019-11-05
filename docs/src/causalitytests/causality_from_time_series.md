
# [Causality from time series](@id causality_time_series)

The `causality` function and its methods provide a common interface for testing causal hypotheses.
For analysing time series, all you need to do is provide a `source` and a `target`. Then, choose 
one of the [available causality tests](@ref causality_tests) to quantify the (directional)
dynamical dependence between `source` and `target`.

For data with uncertainties, see [uncertainty handling](@ref uncertainty_handling).

```@docs
causality(x::AbstractVector, y::AbstractVector, test::CausalityTest)
```
