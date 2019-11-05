
# [Uncertainty handling](@id uncertainty_handling)

All high-level causality test are integrated with the uncertainty handling machinery in 
[UncertainData.jl](https://github.com/kahaaga/UncertainData.jl). Any combination of real-valued vectors, `Vector{<:AbstractUncertainValue}`, or `AbstractUncertainValueDataset` are accepted as inputs to `causality`, making uncertainty quantification on the causality statistics a breeze.

For more fine grained control over the analysis, check out the [syntax overview](@ref syntax_overview) for low-level estimators.

## List of `causality` methods for uncertain data

If both indices and values have uncertainties, use one of the following methods:

- [Naive resampling](@ref causality_uncertain_naiveresampling)
- [Naive constrained resampling](@ref causality_uncertain_naiveconstrained_resampling)
- [Binned resampling](@ref causality_uncertain_binneddatacausalitytest)
- [Strictly increasing, interpolated resampling](@ref causality_uncertain_strictlyincreasing_interpolated)
- [Interpolate-and-bin resampling](@ref causality_uncertain_interpolateandbin_resampling)
