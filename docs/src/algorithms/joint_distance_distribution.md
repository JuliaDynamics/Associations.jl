
## Distribution of joint distances

```@docs
joint_distance_distribution(::Any, ::Any)
```

## Hypothesis test on the joint distance distribution

For the joint distance distribution to indicate a causal influence, it must be significantly 
skewed towards positive values. 

Providing the `OneSampleTTest` type as the first 
argument yields a one sample t-test on the joint distance distribution. From this test, you can extract p-values and obtain 
confidence intervals like in [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) as usual.

```@docs
joint_distance_distribution(::Type{OneSampleTTest}, ::Any, ::Any)
```
