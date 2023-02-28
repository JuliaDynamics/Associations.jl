using StateSpaceSets: AbstractStateSpaceSet, StateSpaceSet
using StateSpaceSets: dimension

export transferentropy
export TransferEntropy
export TransferEntropyEstimator

include("embedding.jl")
include("utils.jl")

"""
    TransferEntropy <: AssociationMeasure

The supertype of all transfer entropy measures. Concrete subtypes are
- [`TEShannon`](@ref)
- [`TERenyiJizba`](@ref)
"""
abstract type TransferEntropy{E, EMB} <: DirectedAssociationMeasure end

"""
The supertype of all dedicated transfer entropy estimators.
"""
abstract type TransferEntropyEstimator end

"""
    transferentropy([measure::TEShannon], est, s, t, [c])
    transferentropy(measure::TERenyiJizba, est, s, t, [c])

Estimate the transfer entropy ``TE^*(S \\to T)`` or ``TE^*(S \\to T | C)`` if `c` is given,
using the provided estimator `est`, where ``*`` indicates the given `measure`.
If `measure` is not given, then `TEShannon(; base = 2)` is the default.

## Arguments

- **`measure`**: The transfer entropy measure, e.g. [`TEShannon`](@ref) or
    [`TERenyi`](@ref), which dictates which formula is computed.
    Embedding parameters are stored in `measure.embedding`, and
    is represented by an [`EmbeddingTE`](@ref) instance. If calling `transferentropy`
    without giving `measure`, then the embedding is optimized by finding
    suitable delay embedding parameters using the ["traditional"](https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/traditional/)
    approach from DynamicalSystems.jl.
- **`s`**: The source timeseries.
- **`t`**: The target timeseries.
- **`c`**: Optional. A conditional timeseries.

## Description

The Shannon transfer entropy is defined as ``TE^S(S \\to T | C) := I^S(T^+; S^- | T^-, C^-)``,
where ``I^S(T^+; S^- | T^-, C^-)`` is [`CMIShannon`](@ref), and marginals for
the CMI are constructed as described in [`EmbeddingTE`](@ref). The definition is
analogous for [`TERenyiJizba`](@ref).

If `s`, `t`, and `c` are univariate timeseries, then the
the marginal embedding variables ``T^+`` (target future), ``T^-`` (target present/past),
``S^-`` (source present/past) and ``C^-`` (present/past of conditioning variables)
are constructed by first jointly embedding  `s`, `t` and `c` with relevant delay
embedding parameters, then subsetting relevant columns of the embedding.

Since estimates of ``TE^*(S \\to T)`` and ``TE^*(S \\to T | C)`` are just a special cases of
conditional mutual information where input data are marginals of a particular form of
[delay embedding](https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/reconstruction/),
*any* combination of variables, e.g. ``S = (A, B)``, ``T = (C, D)``,
``C = (D, E, F)`` are valid inputs (given as `StateSpaceSet`s).
In practice, however, `s`, `t` and `c` are most often timeseries, and if
 `s`, `t` and `c` are [`StateSpaceSet`](@ref)s, it is assumed that the data are
pre-embedded and the embedding step is skipped.

## Compatible estimators


`transferentropy` is just a simple wrapper around [`condmutualinfo`](@ref) that constructs
an appropriate delay embedding from the input data before CMI is estimated. Consequently,
any estimator that can be used for [`ConditionalMutualInformation`](@ref) is, in principle,
also a valid transfer entropy estimator. [`TransferEntropyEstimator`](@ref)s are the
exception - they compute transfer entropy directly.

| Estimator                        | Type                                            | Principle           | [`TEShannon`](@ref) | [`TERenyiJizba`](@ref) |
| -------------------------------- | ----------------------------------------------- | ------------------- | :-----------------: | :--------------------: |
| [`CountOccurrences`](@ref)       | [`ProbabilitiesEstimator`](@ref)                | Frequencies         |         ✓          |           ✓            |
| [`ValueHistogram`](@ref)         | [`ProbabilitiesEstimator`](@ref)                | Binning (histogram) |         ✓          |           ✓            |
| [`Dispersion`](@ref)             | [`ProbabilitiesEstimator`](@ref)                | Dispersion patterns |         ✓          |           ✖            |
| [`Kraskov`](@ref)                | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |           ✖            |
| [`Zhu`](@ref)                    | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |           ✖            |
| [`ZhuSingh`](@ref)               | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |           ✖            |
| [`Gao`](@ref)                    | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |           ✖            |
| [`Goria`](@ref)                  | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |           ✖            |
| [`Lord`](@ref)                   | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |           ✖            |
| [`LeonenkoProzantoSavani`](@ref) | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |           ✓            |
| [`GaussanMI`](@ref)              | [`MutualInformationEstimator`](@ref)            | Parametric          |         ✓          |           ✖            |
| [`KSG1`](@ref)                   | [`MutualInformationEstimator`](@ref)            | Continuous          |         ✓          |           ✖            |
| [`KSG2`](@ref)                   | [`MutualInformationEstimator`](@ref)            | Continuous          |         ✓          |           ✖            |
| [`GaoKannanOhViswanath`](@ref)   | [`MutualInformationEstimator`](@ref)            | Mixed               |         ✓          |           ✖            |
| [`GaoOhViswanath`](@ref)         | [`MutualInformationEstimator`](@ref)            | Continuous          |         ✓          |           ✖            |
| [`FPVP`](@ref)                   | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |         ✓          |           ✖            |
| [`MesnerShalisi`](@ref)          | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |         ✓          |           ✖            |
| [`Rahimzamani`](@ref)            | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |         ✓          |           ✖            |
| [`Zhu1`](@ref)                   | [`TransferEntropyEstimator`](@ref)              | Nearest neighbors   |         ✓          |           ✖            |
| [`Lindner`](@ref)                | [`TransferEntropyEstimator`](@ref)              | Nearest neighbors   |         ✓          |           ✖            |
| [`Hilbert`](@ref)                | [`TransferEntropyEstimator`](@ref)              | Hilbert transform   |         ✓          |           ✖            |
| [`SymbolicTransferEntropy`](@ref)| [`TransferEntropyEstimator`](@ref)              | Hilbert transform   |         ✓          |           ✖            |

"""
function transferentropy end


const TE_ESTIMATORS = Union{
    TransferEntropyEstimator,
    ConditionalMutualInformationEstimator,
    MutualInformationEstimator,
    DifferentialEntropyEstimator,
    ProbabilitiesEstimator,
}

# Embedding optimization
include("optimization/optimization.jl")

include("TEShannon.jl")
include("TERenyiJizba.jl")

function transferentropy(args...; kwargs...)
    return estimate(args...; kwargs...)
end

function estimate(est::TE_ESTIMATORS, args...; kwargs...)
    estimate(TEShannon(), est, args...; kwargs...)
end

function estimate(measure::TransferEntropy, est::TE_ESTIMATORS, x...)
    # If a conditional input (x[3]) is not provided, then C is just a 0-dimensional
    # StateSpaceSet. The horizontal concatenation of C with T then just returns T.
    # We therefore don't need separate methods for the conditional and non-conditional
    # cases.
    S, T, T⁺, C = individual_marginals_te(measure.embedding, x...)
    cmi = te_to_cmi(measure)
    # TE(s -> t) := I(t⁺; s⁻ | t⁻, c⁻).
    return condmutualinfo(cmi, est, T⁺, S, StateSpaceSet(T, C))
end

# When using any estimator except dedicatd `TransferEntropyEstimator`s,
# we use the conditional mutual information decomposition, so we need
# to change the measure for dispatch to work.
te_to_cmi(measure::TEShannon) = CMIShannon(measure.e)
te_to_cmi(measure::TERenyiJizba) = CMIRenyiJizba(measure.e)

function individual_marginals_te(emb::EmbeddingTE, x::AbstractVector...)
    joint, vars, τs, js = te_embed(emb, x...)
    S = joint[:, vars.S]
    T = joint[:, vars.T]
    Tf = joint[:, vars.Tf]
    C = joint[:, vars.C]
    return S, T, Tf, C
end

function h4_marginals(measure::TransferEntropy, x...)
    S, T, T⁺, C = individual_marginals_te(measure.embedding, x...)
    joint = StateSpaceSet(S, T, T⁺, C)
    ST = StateSpaceSet(S, T, C)
    TT⁺ = StateSpaceSet(T, T⁺, C)
    T = StateSpaceSet(T, C)
    return joint, ST, TT⁺, T
end

include("estimators/estimators.jl")
include("convenience/convenience.jl")

# Default to Shannon-type base 2 transfer entropy
function estimate(est::TransferEntropyEstimator, x...)
    estimate(TEShannon(base = 2), est, x...)
end

transferentropy(emb::EmbeddingTE, args...; kwargs...) =
    transferentropy(TEShannon(; embedding = emb), args...; kwargs...)
