# Transfer operator estimation

## Rectangular partitions

```@docs
transferoperator(points, binning_scheme::RectangularBinning; kwargs...)
```

## Triangulated partitions

```@docs
transferoperator(pts, ϵ::TriangulationBinning, simplex_intersection_type::ExactIntersection)
```

```@docs
transferoperator(pts, ϵ::TriangulationBinning, simplex_intersection_type::ApproximateIntersection; n::Int = 200, sample_randomly::Bool = false)
```