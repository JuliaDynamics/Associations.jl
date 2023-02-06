
# [Cross mapping API](@id cross_mapping_api)

This page outlines the cross-mapping API. For concrete implementations of
cross-map based association measures, see [](@ref).

Several cross mapping methods have emerged in the literature
Following Sugihara et al. (2012)'s paper on the convergent cross mapping.
In CausalityTools.jl, we provide a unified interface for using these cross mapping methods.
We indicate the different types of cross mappings by
passing an [`CrossmapMeasure`](@ref) instance as the first argument to [`crossmap`](@ref)
or [`predict`](@ref).

## API

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

### Measures

```@docs
CrossmapMeasure
```

### Estimators

```@docs
CrossmapEstimator
RandomVectors
RandomSegment
ExpandingSegment
```
