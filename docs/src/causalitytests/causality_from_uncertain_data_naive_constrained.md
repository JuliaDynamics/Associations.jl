# [Naive, constrained resampling](@id causality_uncertain_naiveconstrained_resampling)

If you need to truncate the furnishing distributions of your uncertain datasets before 
applying a causality test, use the following method.

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