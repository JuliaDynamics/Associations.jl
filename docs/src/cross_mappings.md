# Cross-mappings

Several cross mapping methods have emerged in the literature
Following Sugihara et al. (2012)'s paper on the convergent cross mapping.
In CausalityTools.jl, we provide a unified interface for using these cross mapping methods.
We indicate the different types of cross mappings by
passing an [`CrossmapMeasure`](@ref) instance as the first argument to [`crossmap`](@ref)
or [`predict`](@ref).

For example, to make predictions `ŷ` for the timeseries `y`, based on an embedding
of timeseries `x`, using the convergent cross mapping method,
do `predict(ConvergentCrossMapping(d = 2), y, x)`.

## Implementations

Any of the estimators below can be used with [`predict`](@ref) and [`crossmap`](@ref).

## API

The cross-mapping API consists of the following two functions:

- [`predict`](@ref), and
- [`crossmap`](@ref).

```@docs
crossmap
predict
CrossmapMeasure
```

## Cross-map measures

```@docs
ConvergentCrossMapping
PairwiseAsymmetricInference
PredictiveDistanceCorrelation
```

## Estimators

```@docs
ExpandingSegment
RandomSegment
```
