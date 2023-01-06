using DelayEmbeddings: AbstractDataset
using ComplexityMeasures: ProbabilitiesEstimator
const VectorOrDataset{D, T} = Union{AbstractVector{T}, AbstractDataset{D, T}} where {D, T}

export CausalityMeasure

"""
    CausalityMeasure

The supertype of all causality measures.
"""
abstract type CausalityMeasure end
