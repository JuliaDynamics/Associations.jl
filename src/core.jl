using DelayEmbeddings: AbstractDataset
using ComplexityMeasures: ProbabilitiesEstimator
const Vector_or_Dataset{D, T} = Union{AbstractVector{T}, AbstractDataset{D, T}} where {D, T}

export CausalityMeasure

"""
    CausalityMeasure

The supertype of all causality measures.
"""
abstract type CausalityMeasure end


# # TODO: all methods shuld use something generic like this.
# function estimate_from_marginals(measure::InformationMeasure, args...;
#         kwargs...) where I <: InformationMeasure
#     msg = "Marginal estimation of $I is not implemented"
#     throw(ArgumentError(msg))
# end

# export estimate
# function estimate end
