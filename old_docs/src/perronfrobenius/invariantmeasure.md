# [Invariant measure estimation](@id invariant_measure_estimation)

## [Rectangular partitions](@id invariant_measure_estimation_rectangular)

```@docs
invariantmeasure(data, binning_scheme::RectangularBinning)
```

## [Triangulated partitions](@id invariant_measure_estimation_triang)

Say we have a 3D delay reconstruction that we have partioned into simplices.

![](triang.png)

There are two methods that approximates invariant measures over that partition.

### Exact simplex intersections

```@docs
invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ExactIntersection)
```

### Approximate simplex intersections

```@docs
invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ApproximateIntersection; n::Int = 100, sample_randomly::Bool = false)
```