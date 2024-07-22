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

## Estimation

- [Example 1](@ref example_TEShannon_EntropyDecomposition_TransferOperator): [`EntropyDecomposition`](@ref) with [`TransferOperator`](@ref) outcome space.
"""
struct TEShannon{B, EMB} <: TransferEntropy
    base::B
    embedding::EMB
    function TEShannon(; base::B = 2, embedding::EMB = EmbeddingTE()) where {B, EMB}
        return new{B, EMB}(base, embedding)
    end
    # TODO: add constructor that automatically determines the embedding.
end

function convert_to_cmi_estimator(est::EntropyDecomposition{<:TEShannon, <:DiscreteInfoEstimator})
    (; definition, est, discretization, pest) = est
    base = definition.base
    return EntropyDecomposition(CMIShannon(; base), est, discretization, pest)
end

function convert_to_cmi_estimator(est::EntropyDecomposition{<:TEShannon, <:DifferentialInfoEstimator})
    return EntropyDecomposition(CMIShannon(; est.definition.base), est.est)
end

function convert_to_cmi_estimator(est::MIDecomposition{<:TEShannon})
    base = est.definition.base
    return MIDecomposition(CMIShannon(; base), est.est)
end

function convert_to_cmi_estimator(est::CMIDecomposition{<:TEShannon})
    base = est.definition.base
    return CMIDecomposition(CMIShannon(; base), est.est)
end


# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
# These are some possible decompositions
# TE(s -> t | c) =
# = I(t⁺; s⁻ | t⁻, c⁻)
# = I(t⁺; s⁻, t⁻, c⁻) - I(t⁺; t⁻, c⁻)
# = h(t⁺ | t⁻,c⁻) - h(t⁺ | s⁻,t⁻,c⁻)
# = h(t⁺, t⁻,c⁻) - h(t⁻,c⁻) - h(t⁺,s⁻,t⁻,c⁻) + h(s⁻,t⁻,c⁻)"

function decomposition_string(
        definition::TEShannon, 
        est::EntropyDecomposition{M, <:DiscreteInfoEstimator{<:Shannon}}
    ) where M
    return "TEₛ(s → t | c) = Hₛ(t⁺, t⁻,c⁻) - Hₛ(t⁻,c⁻) - Hₛ(t⁺,s⁻,t⁻,c⁻) + Hₛ(s⁻,t⁻,c⁻)"
end

function decomposition_string(
    definition::TEShannon, 
    est::EntropyDecomposition{M, <:DifferentialInfoEstimator{<:Shannon}}
    ) where M
    return "TEₛ(s → t | c) = hₛ(t⁺, t⁻,c⁻) - hₛ(t⁻,c⁻) - hₛ(t⁺,s⁻,t⁻,c⁻) + hₛ(s⁻,t⁻,c⁻)"
end

function decomposition_string(
        definition::TEShannon, 
        est::MIDecomposition{M, <:MutualInformationEstimator{<:MIShannon}}
    ) where M
    return "TEₛ(s → t | c) = Iₛ(t⁺; s⁻, t⁻, c⁻) - Iₛ(t⁺; t⁻, c⁻)"
end


function decomposition_string(
        definition::TEShannon, 
        est::CMIDecomposition{M, <:ConditionalMutualInformationEstimator{<:CMIShannon}}
    ) where M
    return "TEₛ(s → t | c) = Iₛ(t⁺; s⁻ | t⁻, c⁻)"
end
