using DelayEmbeddings: AbstractStateSpaceSet
using ComplexityMeasures: ProbabilitiesEstimator
const VectorOrStateSpaceSet{D, T} = Union{AbstractVector{T}, AbstractStateSpaceSet{D, T}} where {D, T}
const ArrayOrStateSpaceSet{D, T, N} = Union{AbstractArray{T, N}, AbstractStateSpaceSet{D, T}} where {D, T, N}

export estimate
export AssociationMeasure
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


# Just use ComplexityMeasures.convert_logunit when it is released.
"""
    _convert_logunit(h_a::Real, , to) → h_b

Convert a number `h_a` computed with logarithms to base `a` to an entropy `h_b` computed
with logarithms to base `b`. This can be used to convert the "unit" of an entropy.
"""
function _convert_logunit(h::Real, base_from, base_to)
    h / log(base_from, base_to)
end
