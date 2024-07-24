"""
    CorrelationMeasure <: AssociationMeasure end

The supertype for correlation measures.

## Concrete implementations

- [`PearsonCorrelation`](@ref)
- [`PartialCorrelation`](@ref)
- [`DistanceCorrelation`](@ref)
"""
abstract type CorrelationMeasure <: AssociationMeasure end

# Future proof definition, to obey the overall API ("estimator contains measure"). 
# Implementations must have `definition` as the first field.
abstract type CorrelationMeasureEstimator{M} <: AssociationMeasure end

include("pearson_correlation.jl")
include("partial_correlation.jl")
include("distance_correlation.jl")
