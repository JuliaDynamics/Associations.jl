# Transfer entropy (TE) estimators]

## [Convenience methods for TE estimation](@id wrapper_TE)

The `transferentropy(driver, response; kwargs...)` and `transferentropy(driver, response, cond; kwargs...)` methods are nice for initial exploration of your data.

*Note: these wrappers use default values for estimator parameters that
may not be suited for your data. For real applications, it is highly
recommended to use the underlying estimators directly. This means that you have to 
perform a [custom delay reconstruction](@ref custom_delay_reconstruction) where 
you explicitly keep track of positions and embeddings delays*, then create a 
[`TEVars`](@ref TEVars) instance that maps the variables of the delay embedding 
to the correct marginals during TE computation.

To understand what is going on under the hood of these methods (strongly recommended), 
see the [general workflow for TE estimation](@ref general_workflow_te).

## Documentation 

### [Regular TE](@id wrapper_regular_TE)

```@docs
transferentropy(::AbstractArray{<:Real, 1}, ::AbstractArray{<:Real, 1})
```

### [Conditional TE](@id wrapper_conditional_TE)

```@docs
transferentropy(::AbstractArray{<:Real, 1}, ::AbstractArray{<:Real, 1}, ::AbstractArray{<:Real, 1})
```