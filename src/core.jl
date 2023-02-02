using DelayEmbeddings: AbstractDataset
using ComplexityMeasures: ProbabilitiesEstimator
const VectorOrDataset{D, T} = Union{AbstractVector{T}, AbstractDataset{D, T}} where {D, T}
const ArrayOrDataset{D, T, N} = Union{AbstractArray{T, N}, AbstractDataset{D, T}} where {D, T, N}

"""
    AssociationMeasure

The supertype of all association measures.
"""
abstract type AssociationMeasure end

# For measures without dedicated estimators, skip the estimator.
function estimate(measure::M, est::Nothing, args...; kwargs...) where M
    estimate(measure, args...; kwargs...)
end

include("contingency/contingency_matrices.jl")
