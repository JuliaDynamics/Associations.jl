using StateSpaceSets: AbstractDataset, Dataset
export transferentropy
export TransferEntropy
export TransferEntropyEstimator

include("embedding.jl")
include("utils.jl")

"""
The supertype of all transfer entropy measures.
"""
abstract type TransferEntropy <: InformationMeasure end

"""
The supertype of all dedicated transfer entropy estimators.
"""
abstract type TransferEntropyEstimator end

"""
    transferentropy(measure::TEShannon, est, s, t, [c])
    transferentropy(measure::TERenyiJizba, est, s, t, [c])

Estimate the transfer entropy ``TE^*(S \\to T)`` or ``TE^*(S \\to T | C)`` if `c` is given,
using the provided estimator `est`, where ``*`` indicates the given `measure`.

## Arguments

- **`measure`**: The transfer entropy measure, e.g. [`TEShannon`](@ref) or
    [`TERenyi`](@ref), which dictates which formula is computed.
    Embedding parameters are stored in `measure.embedding`, and
    is represented by an [`EmbeddingTE`](@ref) instance.
- **`s`**: The source timeseries.
- **`t`**: The target timeseries.
- **`c`**: Optional. Any conditional timeseries.

## Description

The Shannon CMI is defined as ``TE^S(S \\to T | C) &:= I^S(T^+; S^- | T^-, C^-)``,
where ``I^S(T^+; S^- | T^-, C^-)`` is [`CMIShannon`](@ref). The definition is
analogous for [`TERenyiJizba`](@ref).

If `s`, `t`, and `c` are univariate timeseries, then the
the marginal embedding variables ``T^+`` (target future), ``T^-`` (target present/past),
``S^-`` (source present/past) and ``C^-`` (present/past of conditioning variables)
are constructed by first jointly embedding  `s`, `t` and `c` with relevant delay
embedding parameters, then subsetting relevant columns of the embedding.
Details are found in [`EmbeddingTE`](@ref).

Since estimates of ``TE^*(S \\to T)`` and ``TE^*(S \\to T | C)`` are just a special cases of
conditional mutual information where input data are marginals of a particular form of
[delay embedding](https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/reconstruction/),
*any* combination of variables, e.g. ``S = (A, B)``, ``T = (C, D)``,
``C = (D, E, F)`` are valid inputs (given as `Dataset`s).
In practice, however, `s`, `t` and `c` are most often timeseries, and if
 `s`, `t` and `c` are [`Dataset`](@ref)s, it is assumed that the data are
pre-embedded and the embedding step is skipped.

## Compatible estimators

`transferentropy` is just a simple wrapper around [`condmutualinfo`](@ref) that constructs
an appropriate delay embedding from the input data before CMI is estimated. Consequently,
any estimator that can be used for [`ConditionalMutualInformation`](@ref) is, in principle,
also a valid transfer entropy estimator. Documentation strings for [`TEShannon`](@ref) and
[`TERenyiJizba`](@ref)) list compatible estimators, and an overview table can be found in
the online documentation.
"""
function transferentropy end

include("TEShannon.jl")
include("TERenyiJizba.jl")

function transferentropy(measure::TransferEntropy, est, x...)
    # If a conditional input (x[3]) is not provided, then C is just a 0-dimensional
    # dataset. The horizontal concatenation of C with T then just returns T.
    # We therefore don't need separate methods for the conditional and non-conditional
    # cases.
    S, T, 𝒯, C = individual_marginals(measure.embedding, x...)
    cmi = convert_to_cmi_measure(measure)
    # TE(s -> t) := I(t⁺; s⁻ | t⁻, c⁻).
    return condmutualinfo(cmi, est, 𝒯, S, Dataset(T, C))
end

convert_to_cmi_measure(measure::TEShannon) = CMIShannon(measure.e)
convert_to_cmi_measure(measure::TERenyiJizba) = CMIRenyiJizba(measure.e)

function individual_marginals(emb::EmbeddingTE, x::AbstractVector...)
    joint, vars, τs, js = te_embed(emb, x...)
    S = joint[:, vars.S]
    T = joint[:, vars.T]
    𝒯 = joint[:, vars.𝒯]
    C = joint[:, vars.C]
    return S, T, 𝒯, C
end

function h4_marginals(measure::TransferEntropy, x...)
    S, T, T⁺, C = individual_marginals(measure.embedding, x...)
    joint = Dataset(S, T, T⁺, C)
    ST = Dataset(S, T, C)
    TT⁺ = Dataset(T, T⁺, C)
    T = Dataset(T, C)
    return joint, ST, TT⁺, T
end

include("estimators/estimators.jl")
