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

function association(est::JointProbabilities{<:SECMI}, x, y, z::Vector{AbstractVector}...)
    probs = probabilities(est.discretization, x, y, z...)
    return association(est.definition, probs)
end

function association(est::JointProbabilities{<:SECMI}, x, y, z::AbstractStateSpaceSet)
    probs = probabilities(est.discretization, x, y, columns(z)...)
    return association(est.definition, probs)
end

"""
    SECMITest <: IndependenceTest
    SECMITest(; nshuffles = 19, surrogate = RandomShuffle(), rng = Random.default_rng())

A test for conditional independence based on the [`SECMI`](@ref) measure.
"""
Base.@kwdef struct SECMITest{E, S, I, RNG} <: IndependenceTest{E}
    est_or_measure::E = SECMI(; base = 2)
    surrogate::S = RandomShuffle()
    nshuffles::I = 19
    rng::RNG = Random.default_rng()
end

function independence(test::SECMITest, x, y, z)
    (; est_or_measure, surrogate, nshuffles, rng) = test
    sx = surrogenerator(x, surrogate, rng)
    secmiₖ = zeros(test.n)
    for k = 1:nshuffles
        secmiₖ[k] = association(est_or_measure, sx(), y, z)
    end
    μ̂ = 1/nshuffles * sum(secmiₖ)
    σ̂ = 1/(nshuffles - 1) * sum((sₖ - μ̂)^2 for sₖ in secmiₖ)
    emp_cdf::Function = ecdf(secmiₖ) 
    𝒩 = Normal(μ̂, σ̂)
    # degrees of freedom for Chi squared distribution estimated as the mean of the `secmiₖ`.
    # (page 18 in Kubkowski et al)
    𝒳 = Chisq(μ̂)

end