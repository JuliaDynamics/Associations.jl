using TimeseriesSurrogates
using Random
using StatsBase: ecdf
using Distributions

export ShortExpansionConditionalMutualInformation, SECMI

"""
    ShortExpansionConditionalMutualInformation <: MultivariateInformationMeasure
    ShortExpansionConditionalMutualInformation(; base = 2)
"""
Base.@kwdef struct ShortExpansionConditionalMutualInformation{B} <: MultivariateInformationMeasure
    base::B = 2
end
const SECMI = ShortExpansionConditionalMutualInformation

function Base.show(io::IO, m::ShortExpansionConditionalMutualInformation)
    msg = "ShortExpansionConditionalMutualInformation(; base = $(m.base))"
    print(io, msg)
end

# Assumes 1st dimension of `probs` corresponds to X, 2nd dimension of `probs`
# corresponds to Y, and dimensions `3:ndims(probs)` correspond to marginals Zₖ, 
function association(definition::SECMI, probs::Probabilities{T, N}) where {T, N}
    @assert N >= 3
    (; base) = definition

    def_mi = MIShannon(; base = base)
    def_cmi = CMIShannon(; base = base)

    m = ndims(probs) - 2
    pXY = marginal(probs, dims = 1:2)
    mi_XY = association(def_mi, pXY)
    cmis = 0.0
    for k = 1:m
        dims = (1, 2, 2 + k)
        cmis += association(def_cmi, marginal(probs, dims = dims))
    end

    return (1 - m) * mi_XY + cmis
end

function association(est::JointProbabilities{<:SECMI}, x, y, z...)
    probs = probabilities(est.discretization, x, y, z...)
    return association(est.definition, probs)
end

function association(est::JointProbabilities{<:SECMI}, x, y, z::AbstractStateSpaceSet)
    probs = probabilities(est.discretization, x, y, columns(z)...)
    return association(est.definition, probs)
end
