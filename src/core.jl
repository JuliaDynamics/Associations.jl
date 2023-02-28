using DelayEmbeddings: AbstractStateSpaceSet
using ComplexityMeasures: ProbabilitiesEstimator
const VectorOrStateSpaceSet{D, T} = Union{AbstractVector{T}, AbstractStateSpaceSet{D, T}} where {D, T}
const ArrayOrStateSpaceSet{D, T, N} = Union{AbstractArray{T, N}, AbstractStateSpaceSet{D, T}} where {D, T, N}

export DirectedAssociationMeasure
"""
    AssociationMeasure

The supertype of all association measures.
"""
abstract type AssociationMeasure end

abstract type DirectedAssociationMeasure <: AssociationMeasure end

# For measures without dedicated estimators, skip the estimator.
function estimate(measure::M, est::Nothing, args...; kwargs...) where M
    estimate(measure, args...; kwargs...)
end

include("contingency_matrices.jl")
