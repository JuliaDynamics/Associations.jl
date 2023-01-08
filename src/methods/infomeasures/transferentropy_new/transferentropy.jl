using StateSpaceSets: AbstractDataset, Dataset
export transferentropy

include("utils.jl")

"""
The supertype of all transfer entropy measures.
"""
abstract type TransferEntropy <: InformationMeasure end

"""
    transferentropy(measure::TEShannon, est, s, t, c; kwargs...)
    transferentropy(measure::TERenyiJizba, est, s, t, c; kwargs...)

Estimate the given transfer entropy `measure` from source variable ``S`` to
target variable ``T``, conditioned on conditional
variable(s) ``C`` using the given estimator `est`.

This is just a simple wrapper around [`condmutualinfo`](@ref), and `est` can be
any [`ConditionalMutualInformationEstimator`](@ref), [`MutualInformationEstimator`](@ref),
,[`DifferentialEntropyEstimator`](@ref), or any [`ProbabilitiesEstimator`](@ref) that
accepts multivariate input data or has an implementation for [`marginal_encodings`](@ref).

## Description

Transfer entropy is essentially just a special case of conditional mutual information
where the input data is a certain type of delay embedding. Here, we use
the [`TEShannon`](@ref) measure to illustrate the embedding procedure,
which is also used for [`TERenyiJizba`](@ref).

The Shannon CMI is defined as ``TE^S(S \\to T | C) &:= I^S(T^+; S^- | T^-, C^-)``,
where the variables ``T^+`` (target future), ``T^-``
(target present/past), ``S^-`` (source present/past) and ``C^-`` (present/past
of conditioning variables) are constructed by first jointly embedding
``S``, ``T`` and ``C`` with relevant delay embedding parameters, then subsetting
relevant columns of the embedding.

Since ``TE^S(S \\to T)`` is just a special case of conditional mutual information,
*any* combination of variables, e.g. ``S = (A, B)``, ``T = (C, D)``, ``C = (D, E, F)`` are
valid inputs. In practice, however, the curse of dimensionality quickly
slows down computation and reliability of the estimated transfer entropy,
so typically the input data are *timeseries*.

"""
function transferentropy end

include("TEShannon.jl")
include("TERenyiJizba.jl")

function transferentropy(measure::TransferEntropy, est, x...;
        τT = -1, τS = -1, η𝒯 = 1, dT = 1, dS = 1, d𝒯 = 1,τC = -1, dC = 1)
    # If a conditional input (x[3]) is not provided, then C is just a 0-dimensional
    # dataset. The horizontal concatenation of C with T then just returns T.
    # We therefore don't need separate methods for the conditional and non-conditional
    # cases.
    emb = EmbeddingTE(; τT, τS, τC, η𝒯, dT, dS, d𝒯, dC)
    S, T, 𝒯, C = te_marginals(emb, x...)
    cmi = get_cmi_measure(measure)
    # TE(s -> t) := I(t⁺; s⁻ | t⁻, c⁻).
    return condmutualinfo(cmi, est, 𝒯, S, Dataset(T, C))
end

get_cmi_measure(measure::TEShannon) = CMIShannon(measure.e)
get_cmi_measure(measure::TERenyiJizba) = CMIRenyiJizba(measure.e)

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
