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
"""
struct TEShannon{B, EMB} <: TransferEntropy{B, EMB}
    base::B
    embedding::EMB
    function TEShannon(; base::B = 2, embedding::EMB = EmbeddingTE()) where {B, EMB}
        return new{B, EMB}(base, embedding)
    end
    # TODO: add constructor that automatically determines the embedding.
end

function convert_to_cmi_estimator(est::MIDecomposition{<:TEShannon})
    base = est.definition.base
    return MIDecomposition(CMIShannon(; base), est.est)
end

function convert_to_cmi_estimator(est::EntropyDecomposition{<:TEShannon})
    base = est.definition.base
    return EntropyDecomposition(CMIShannon(; base), est.est)
end

function convert_to_cmi_estimator(est::CMIDecomposition{<:TEShannon})
    base = est.definition.base
    return CMIDecomposition(CMIShannon(; base), est.est)
end
