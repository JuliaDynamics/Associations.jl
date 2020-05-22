# [Interpolate-and-bin resampling](@id causality_uncertain_interpolateandbin_resampling)

!!! note
    Uncertain data handling relies on the [UncertainData](https://github.com/kahaaga/UncertainData.jl).
    To use the resampling schemes, you need to load this package by first running `using UncertainData` 
    in the Julia console.

## InterpolateBinTest

```@docs
InterpolateBinTest
```

## Call signature

```@docs
causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::InterpolateBinTest)
```
