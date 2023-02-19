using DelayEmbeddings: delay_f1nn
export TEShannon

"""
    TEShannon <: TransferEntropy
    TEShannon(; base = 2; embedding = EmbeddingTE()) <: TransferEntropy

The Shannon-type transfer entropy measure.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise
    and conditional dependence.
- Use with [`transferentropy`](@ref) to compute the raw transfer entropy.

## Description

The transfer entropy from source ``S`` to target ``T``, potentially
conditioned on ``C`` is defined as

```math
\\begin{align*}
TE(S \\to T) &:= I^S(T^+; S^- | T^-) \\\\
TE(S \\to T | C) &:= I^S(T^+; S^- | T^-, C^-)
\\end{align*}
```

where ``I(T^+; S^- | T^-)`` is the Shannon conditional mutual information
([`CMIShannon`](@ref)). The variables ``T^+``, ``T^-``,
``S^-`` and ``C^-`` are described in the docstring for [`transferentropy`](@ref).

## Compatible estimators

Shannon-type transfer entropy can be estimated using a range of different estimators,
which all boil down to computing conditional mutual information, except for
[`TransferEntropyEstimator`](@ref), which compute transfer entropy using some direct method.

| Estimator                        | Type                                            | Principle           | [`TEShannon`](@ref) |
| -------------------------------- | ----------------------------------------------- | ------------------- | :-----------------: |
| [`CountOccurrences`](@ref)       | [`ProbabilitiesEstimator`](@ref)                | Frequencies         |         ✓          |
| [`ValueHistogram`](@ref)         | [`ProbabilitiesEstimator`](@ref)                | Binning (histogram) |         ✓          |
| [`SymbolicPermuation`](@ref)     | [`ProbabilitiesEstimator`](@ref)                | Ordinal patterns    |         ✓          |
| [`Dispersion`](@ref)             | [`ProbabilitiesEstimator`](@ref)                | Dispersion patterns |         ✓          |
| [`Kraskov`](@ref)                | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |
| [`Zhu`](@ref)                    | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |
| [`ZhuSingh`](@ref)               | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |
| [`Gao`](@ref)                    | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |
| [`Goria`](@ref)                  | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |
| [`Lord`](@ref)                   | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |
| [`LeonenkoProzantoSavani`](@ref) | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |         ✓          |
| [`GaussanMI`](@ref)              | [`MutualInformationEstimator`](@ref)            | Parametric          |         ✓          |
| [`KSG1`](@ref)                   | [`MutualInformationEstimator`](@ref)            | Continuous          |         ✓          |
| [`KSG2`](@ref)                   | [`MutualInformationEstimator`](@ref)            | Continuous          |         ✓          |
| [`GaoKannanOhViswanath`](@ref)   | [`MutualInformationEstimator`](@ref)            | Mixed               |         ✓          |
| [`GaoOhViswanath`](@ref)         | [`MutualInformationEstimator`](@ref)            | Continuous          |         ✓          |
| [`FPVP`](@ref)                   | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |         ✓          |
| [`MesnerShalisi`](@ref)          | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |         ✓          |
| [`Rahimzamani`](@ref)            | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |         ✓          |
| [`Zhu1`](@ref)                   | [`TransferEntropyEstimator`](@ref)              | Nearest neighbors   |         ✓          |
| [`Lindner`](@ref)                | [`TransferEntropyEstimator`](@ref)              | Nearest neighbors   |         ✓          |
"""
struct TEShannon{E <: Shannon, EMB} <: TransferEntropy{E, EMB}
    e::E
    embedding::EMB
    function TEShannon(; base = 2, embedding::EMB = EmbeddingTE()) where EMB
        e = Shannon(; base = base)
        return new{typeof(e), EMB}(e, embedding)
    end
    function TEShannon(e::E; embedding::EMB = EmbeddingTE()) where {E <:Shannon, EMB}
        return new{E, EMB}(e, embedding)
    end
    # TODO: add constructor that automatically determines the embedding.
end

# function transferentropy(
#         est::Union{
#             ConditionalMutualInformationEstimator,
#             MutualInformationEstimator,
#             DifferentialEntropyEstimator,
#             ProbabilitiesEstimator
#         },
#         x...; kwargs...)
#     N = length(first(x))

#     # A very naive heuristic to avoid too high dimensions. *All* marginals are optimised,
#     # so in the worst case, the dimension triples.
#     maxdim = floor(Int, N^(1/7))
#     # The maxlag should also scale with the length the input.
#     maxlag = min(floor(Int, N ÷ 50), 100)
#     dmethod = "mi_min"
#     method = delay_f1nn
#     opt = OptimiseTraditional(; maxdim, maxlag, method, dmethod)
#     m = TEShannon(; base = 2, embedding = EmbeddingTE(opt, x...))
#     return transferentropy(m, est, x...; kwargs...)
# end

# If a pre-computed [`ContingencyMatrix`](@ref) `c` is provided, then we just compute
# the conditional mutual information directly from it, assuming the contingency matrix
# was constructed from a meaningful embedding.
function transferentropy(measure::TEShannon, c::ContingencyMatrix)
    cmi = CMIShannon(; base = measure.base)
    return condmutualinfo(cmi, c)
end
