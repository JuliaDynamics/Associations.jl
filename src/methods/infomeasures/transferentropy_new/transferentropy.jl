using StateSpaceSets: AbstractDataset, Dataset
export TEShannon
export transferentropy

include("utils.jl")

abstract type TransferEntropy <: InformationMeasure end

"""
    TEShannon(; base = 2)

The Shannon-type transfer entropy measure given by
``TE(s \\to t | c) := I(t⁺; s⁻ | t⁻, c⁻)``. Used
with [`transferentropy`](@ref).
"""
Base.@kwdef struct TEShannon{B} <: TransferEntropy
    base::B = 2
end

"""
    transferentropy([measure::TEShannon], est, s, t, c; kwargs...)
    transferentropy([measure::TEShannon], c::ContingencyMatrix; kwargs...)

Estimate the transfer entropy ``TE(s \\to t | c) := I(t⁺; s⁻ | t⁻, c⁻)`` for
source variable `s`, target variable `t` and conditional variable(s) `c`.

This is just a simple wrapper around [`condmutualinfo`](@ref), and `est` can be
any [`ConditionalMutualInformationEstimator`](@ref), [`MutualInformationEstimator`](@ref),
,[`DifferentialEntropyEstimator`](@ref), or any [`ProbabilitiesEstimator`](@ref) that
accepts multivariate input data or has an implementation for [`marginal_encodings`](@ref).

If a pre-computed [`ContingencyMatrix`](@ref) `c` is provided, then we just compute
the conditional mutual information directly from it, assuming the contingency matrix
was constructed from a meaningful embedding.
"""
function transferentropy(measure::TEShannon, est, x...;
        τT = -1, τS = -1, η𝒯 = 1, dT = 1, dS = 1, d𝒯 = 1,τC = -1, dC = 1)
    # TE(s -> t) := I(t⁺; s⁻ | t⁻, c⁻).
    # If a conditional input (x[3]) is not provided, then C is just a 0-dimensional
    # dataset. The horizontal concatenation of C with T then just returns T.
    # We therefore don't need separate methods for the conditional and non-conditional
    # cases.
    emb = EmbeddingTE(; τT, τS, τC, η𝒯, dT, dS, d𝒯, dC)
    S, T, 𝒯, C = te_marginals(emb, x...)
    cmi = CMIShannon(; base = measure.base)
    return condmutualinfo(cmi, est, 𝒯, S, Dataset(T, C))
end

function transferentropy(
        est::Union{
            ConditionalMutualInformationEstimator,
            MutualInformationEstimator,
            DifferentialEntropyEstimator,
            ProbabilitiesEstimator
        },
        args...; kwargs...)
    return transferentropy(TEShannon(; base = 2), est, args...; kwargs...)
end

function transferentropy(measure::TEShannon, c::ContingencyMatrix)
    cmi = CMIShannon(; base = measure.base)
    return condmutualinfo(cmi, c)
end

function te_marginals(emb::EmbeddingTE, x::AbstractVector...)
    pts, vars, τs, js = te_embed(emb, x...)
    return te_marginals(vars, pts)
end

# map a set of pre-embedded points to the correct marginals for transfer entropy computation
function te_marginals(vars::TEVars, pts::AbstractDataset)
    S = pts[:, vars.S]
    T = pts[:, vars.T]
    𝒯 = pts[:, vars.𝒯]
    C = pts[:, vars.C]
    return S, T, 𝒯, C
end


#include("hilbert.jl")
include("TERenyi.jl")
