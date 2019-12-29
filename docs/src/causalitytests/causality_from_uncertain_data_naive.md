
# [Testing for causality from uncertain time series](@id causality_uncertain_naiveresampling)

!!! note

    The use the causality tests with uncertain data, first load the 
    [UncertainData](https://github.com/kahaaga/UncertainData.jl) package by running 
    `using UncertainData` in the Julia console.

## Syntax

### Uncertainties only in values

```@docs
causality(source::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset},
        target::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset},
        test::CausalityTest)
```

### Uncertainties in both indices and values

```@docs
causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest)
```
