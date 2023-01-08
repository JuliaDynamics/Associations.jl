export TEShannon

"""
    TEShannon(; base = 2) <: TransferEntropy

The Shannon-type transfer entropy measure.

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

## Estimation

Shannon transfer entropy can be estimated using one of four estimator types

- **[`ProbabilitiesEstimator`](@ref)**. Decomposes the transfer entropy into a sum of
    marginal entropies, then explicitly computes a probability distribution for each
    marginal based
    on some chosen property of the data. Marginal entropies of type `e` are then computed by
    giving these probabilities to [`entropy`](@ref). Example: [`ValueHistogram`](@ref).
- **[`DifferentialDifferentialEntropyEstimator`](@ref)**. Decomposes the transfer entropy
    into a sum of marginal entropies, then computes the differential entropy of type `e`
    for each marginal. Example: [`Kraskov`](@ref).
- **[`MutualInformationEstimator`](@ref)**. Decomposes the transfer entropy into a sum of
    mutual information terms, which are estimated separately using [`mutualinfo`](@ref).
    Example: [`KraskovStögbauerGrassberger2`](@ref).
- **[`ConditionalMutualInformationEstimator`](@ref)**. Use a dedicated CMI estimator.
    Example: [`FrenzelPompeVelmejkaPalus`](@ref).
"""
struct TEShannon{E <: Shannon, EMB} <: TransferEntropy
    e::E
    embedding::EMB
    function TEShannon(; base = 2, embedding::EMB = EmbeddingTE()) where EMB
        e = Shannon(; base = base)
        return new{typeof(e), EMB}(e, embedding)
    end
    function TEShannon(e::E; embedding::EMB = EmbeddingTE()) where {E <:Shannon, EMB}
        return new{E, EMB}(e)
    end
    # TODO: add constructor that automatically determines the embedding.
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

# If a pre-computed [`ContingencyMatrix`](@ref) `c` is provided, then we just compute
# the conditional mutual information directly from it, assuming the contingency matrix
# was constructed from a meaningful embedding.
function transferentropy(measure::TEShannon, c::ContingencyMatrix)
    cmi = CMIShannon(; base = measure.base)
    return condmutualinfo(cmi, c)
end
