
# [Cross mapping API](@id cross_mapping_api)

The cross mapping API define different ways of quantifying association based on the 
concept of "cross mapping", which has appeared in many contexts in the literature,
and gained huge popularity with  [Sugihara2012](@citet)'s on *convergent cross mapping*.

Since their paper, several cross mapping methods and frameworks have emerged in the
literature. In CausalityTools.jl, we provide a unified interface for using these cross
mapping methods.

To estimate a cross map measure, you simply input a [`CrossmapMeasure`](@ref) instance
as the first argument to a [`CrossmapEstimator`](@ref), which is then fed to 
the [`crossmap`](@ref) or [`predict`](@ref) functions. 

The cross mapping API consists of the following functions.

- [`predict`](@ref)
- [`crossmap`](@ref)

These functions can dispatch on a [`CrossmapMeasure`](@ref), and we currently implement

- [`ConvergentCrossMapping`](@ref).
- [`PairwiseAsymmetricEmbedding`](@ref).

```@docs
crossmap
predict
```

## Measures

```@docs
CrossmapMeasure
ConvergentCrossMapping
PairwiseAsymmetricInference
```

## Estimators

```@docs
CrossmapEstimator
RandomVectors
RandomSegment
ExpandingSegment
```
