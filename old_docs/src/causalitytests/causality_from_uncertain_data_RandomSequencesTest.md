
# [Random chunk resampling](@id causality_uncertain_randomsequencestest)

!!! note
    Uncertain data handling relies on the [UncertainData](https://github.com/kahaaga/UncertainData.jl).
    To use the resampling schemes, you need to load this package by first running `using UncertainData` 
    in the Julia console.
    
## RandomSequencesTest

```@docs
RandomSequencesTest
```

## Call signature

```@docs
causality(
        x::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset, AbstractUncertainIndexValueDataset}, 
        y::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset, AbstractUncertainIndexValueDataset},
        test::RandomSequencesTest)
```
