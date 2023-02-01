using DelayEmbeddings: AbstractDataset
using ComplexityMeasures: ProbabilitiesEstimator
const VectorOrDataset{D, T} = Union{AbstractVector{T}, AbstractDataset{D, T}} where {D, T}
const ArrayOrDataset{D, T, N} = Union{AbstractArray{T, N}, AbstractDataset{D, T}} where {D, T, N}

export CausalityMeasure

"""
    CausalityMeasure

The supertype of all causality measures.
"""
abstract type CausalityMeasure end

# For measures without dedicated estimators, skip the estimator.
function estimate(measure::M, est::Nothing, args...; kwargs...) where M
    estimate(measure, args...; kwargs...)
end

include("contingency/contingency_matrices.jl")
