# Invariant measure estimation

## Rectangular partitions

```@docs
invariantmeasure(data, binning_scheme::RectangularBinning)
```

## Triangulated partitions

```@docs
invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ExactIntersection)
```

```@docs
invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ApproximateIntersection; n::Int = 100, sample_randomly::Bool = false)
```