export TERenyiJizba

"""
    TERenyiJizba() <: TransferEntropy

The Rényi transfer entropy from Jizba et al. (2012).

## Description

The transfer entropy from source ``S`` to target ``T``, potentially
conditioned on ``C`` is defined as

```math
\\begin{align*}
TE(S \\to T) &:= I_q^{R_J}(T^+; S^- | T^-) \\\\
TE(S \\to T | C) &:= I_q^{R_J}(T^+; S^- | T^-, C^-),
\\end{align*},
```
where ``I_q^{R_J}(T^+; S^- | T^-)`` is Jizba et al. (2012)'s definition of
conditional mutual information ([`CMIRenyiJizba`](@ref)).
The variables ``T^+``, ``T^-``,
``S^-`` and ``C^-`` are described in the docstring for [`transferentropy`](@ref).

## Estimation

Rènyi transfer entropy of this type can be estimated using one the following estimator
types.

- **[`ProbabilitiesEstimator`](@ref)**. Decomposes the transfer entropy into a sum of
    marginal entropies, then explicitly computes a probability distribution for each
    marginal based
    on some chosen property of the data. Marginal entropies of type `e` are then computed by
    giving these probabilities to [`entropy`](@ref). Example: [`ValueHistogram`](@ref).
- **[`DifferentialDifferentialEntropyEstimator`](@ref)**. Decomposes the transfer entropy
    into a sum of marginal entropies, then computes the differential entropy of type `e`
    for each marginal. Example: [`Kraskov`](@ref).
"""
struct TERenyiJizba{E <: Renyi, EMB} <: TransferEntropy
    e::E
    embedding::EMB
    function TERenyiJizba(; base = 2, q = 1.5, embedding::EMB = EmbeddingTE()) where EMB
        e = Renyi(; base = base, q = q)
        return new{typeof(e), EMB}(e, embedding)
    end
    function TERenyiJizba(e::E; embedding::EMB = EmbeddingTE()) where {E <: Renyi, EMB}
        return new{E, EMB}(e, embedding)
    end
end

"""
    escort_distribution(probs, i::Int, q::Real)

The escort distribution for a probability distribution `probs`. For `q > 1`, the
escort distribution emphasises more probable events and de-emphasises more improbable
events. For `q < 1`, the situation is reversed.

```math
\\text{esc}_q(x) = \\dfrac{p^q(x)}{\\sum_{x \\in \\mathcal{X}} p^q(x)}
```
"""
function escort_distribution(probs, i::Int, q)
    return probs[i]^q / sum(probs .^ q)
end
