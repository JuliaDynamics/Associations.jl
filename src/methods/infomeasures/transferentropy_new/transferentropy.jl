using StateSpaceSets: AbstractDataset, Dataset
export TEShannon
export transferentropy

include("utils.jl")


"""
    TEShannon(; base = 2)

The Shannon-type transfer entropy measure given by
``TE(s \\to t | c) := I(t⁺; s⁻ | t⁻, c⁻)``. Used
with [`transferentropy`](@ref).
"""
Base.@kwdef struct TEShannon{B} <: InformationMeasure
    base::B = 2
end

"""
    transferentropy([measure::TEShannon, est, s, t, c; kwargs...)
    transferentropy([measure::TEShannon, c::ContingencyMatrix; kwargs...)

Estimate the transfer entropy ``TE(s \\to t | c) := I(t⁺; s⁻ | t⁻, c⁻)`` for
source variable `s`, target variable `t` and conditional variable(s) `c`.

This is just a simple wrapper around [`condmutualinfo`](@ref), and `est` can be
any [`ConditionalMutualInformationEstimator`](@ref), [`MutualInformationEstimator`](@ref),
,[`DifferentialEntropyEstimator`](@ref), or any [`ProbabilitiesEstimator`](@ref) that
accepts multivariate input data or has an implementation for [`marginal_encodings`](@ref).
If a pre-computed [`ContingencyMatrix`](@ref) `c` is provided, then compute
transfer entropy directly from it.

"""
function transferentropy(measure::TEShannon, est, args...; kwargs...)
    # TE(s -> t) := I(t⁺; s⁻ | t⁻, c⁻).
    # C may be a 0-dimensional dataset, but that works fine with concatenation etc, so
    # we can use just one common method here.
    S, T, 𝒯, C = te_marginals(EmbeddingTE(; kwargs...), args...)
    cmi = CMIShannon(; base = measure.base)
    return condmutualinfo(cmi, est, 𝒯, S, Dataset(T, C))
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
