
# [Naive resampling](@id causality_uncertain_naiveresampling)

## Uncertainties only in values

```@docs
causality(source::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset},
        target::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset},
        test::CausalityTest)
```

## Uncertainties in both indices and values

```@docs
causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest)
```
