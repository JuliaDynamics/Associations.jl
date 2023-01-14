export TEShannon

"""
    TEShannon <: TransferEntropy
    TEShannon(; base = 2; embedding = EmbeddingTE()) <: TransferEntropy

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

## Compatible estimators

- **[`ProbabilitiesEstimator`](@ref)**: Any probabilities estimator that accepts
    multivariate input data or has an implementation for [`marginal_encodings`](@ref).
    Transfer entropy is computed a sum of marginal (discrete) entropy estimates.
    Example: [`ValueHistogram`](@ref).
- **[`DifferentialEntropyEstimator`](@ref)**. Any differential entropy
    estimator that accepts multivariate input data.
    Transfer entropy is computed a sum of marginal differential entropy estimates.
    Example: [`Kraskov`](@ref).
- **[`MutualInformationEstimator`](@ref)**. Any mutual information estimator.
    Formulates the transfer entropy as a sum of mutual information terms, which are
    estimated separately using [`mutualinfo`](@ref). Example: [`KraskovStögbauerGrassberger2`](@ref).
- **[`ConditionalMutualInformationEstimator`](@ref)**. Dedicated CMI estimators.
    Example: [`FrenzelPompeVelmejkaPalus`](@ref).
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

function transferentropy(
        est::Union{
            ConditionalMutualInformationEstimator,
            MutualInformationEstimator,
            DifferentialEntropyEstimator,
            ProbabilitiesEstimator
        },
        x...; kwargs...)
    emb = EmbeddingTE(OptimiseTraditional(maxdim = 5, maxlag = 0.05), x...)
    m = TEShannon(; base = 2, embedding = emb)
    return transferentropy(m, est, x...; kwargs...)
end

# If a pre-computed [`ContingencyMatrix`](@ref) `c` is provided, then we just compute
# the conditional mutual information directly from it, assuming the contingency matrix
# was constructed from a meaningful embedding.
function transferentropy(measure::TEShannon, c::ContingencyMatrix)
    cmi = CMIShannon(; base = measure.base)
    return condmutualinfo(cmi, c)
end
