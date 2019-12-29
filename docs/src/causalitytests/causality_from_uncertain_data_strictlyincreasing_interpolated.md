# [Strictly increasing, interpolated resampling](@id causality_uncertain_strictlyincreasing_interpolated)

!!! note
    Uncertain data handling relies on the [UncertainData](https://github.com/kahaaga/UncertainData.jl).
    To use the resampling schemes, you need to load this package by first running `using UncertainData` 
    in the Julia console.
    
## Call signature

```@docs
causality(source::AbstractUncertainIndexValueDataset, 
        target::AbstractUncertainIndexValueDataset, 
        test::CausalityTest, constraint::SequentialSamplingConstraint, grid::RegularGrid)
```
