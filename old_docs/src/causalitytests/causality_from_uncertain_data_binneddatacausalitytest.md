
# [Binned resampling](@id causality_uncertain_binneddatacausalitytest)

!!! note
    Uncertain data handling relies on the [UncertainData](https://github.com/kahaaga/UncertainData.jl).
    To use the resampling schemes, you need to load this package by first running `using UncertainData` 
    in the Julia console.

## BinnedDataCausalityTest

```@docs
BinnedDataCausalityTest
```

## Call signature

```@docs
causality(::AbstractUncertainIndexValueDataset, y::AbstractUncertainIndexValueDataset, test::BinnedDataCausalityTest)
```
