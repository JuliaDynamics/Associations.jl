export ClosenessMeasure
"""
    ClosenessMeasure <: AssociationMeasure

The supertype for all multivariate information-based measure definitions.

## Implementations

- [`JointDistanceDistribution`](@ref)
- [`SMeasure`](@ref)
- [`HMeasure`](@ref)
- [`MMeasure`](@ref)
- [`LMeasure`](@ref)

"""
abstract type ClosenessMeasure <: AssociationMeasure end

include("JointDistanceDistribution.jl")
include("SMeasure.jl")
include("HMeasure.jl")
include("MMeasure.jl")
include("LMeasure.jl")
include("common.jl")
