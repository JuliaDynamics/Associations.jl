# [Naive, constrained resampling](@id causality_uncertain_naiveconstrained_resampling)

If you need to truncate the furnishing distributions of your uncertain datasets before 
applying a causality test, use the following method.

!!! note
    Uncertain data handling relies on the [UncertainData](https://github.com/kahaaga/UncertainData.jl).
    To use the resampling schemes, you need to load this package by first running `using UncertainData` 
    in the Julia console.
    

## ConstrainedTest

```@docs
ConstrainedTest
```

## Call signature

```@docs
causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::ConstrainedTest)
```