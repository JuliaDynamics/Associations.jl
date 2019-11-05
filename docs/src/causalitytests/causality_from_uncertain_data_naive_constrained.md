# [Naive, constrained resampling](@id causality_uncertain_naiveconstrained_resampling)

If you need to truncate the furnishing distributions of your uncertain datasets before 
applying a causality test, use the following method.

```@docs
causality(source, target, test::CausalityTest, resampling::ConstrainedResampling)
```

## ConstrainedResampling

```@docs
ConstrainedResampling
```
