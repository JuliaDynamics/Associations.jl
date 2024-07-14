using ComplexityMeasures: Tsallis

export ConditionalEntropyTsallisAbe

"""
    ConditionalEntropyTsallisAbe <: ConditionalEntropy
    ConditionalEntropyTsallisAbe(; base = 2, q = 1.5)

[Abe2001](@citet)'s discrete Tsallis conditional entropy measure.

## Usage 

- Use with [`association`](@ref) to compute the Tsallis-Abe conditional entropy between 
    two variables.

## Compatible estimators

- [`JointProbabilities`](@ref)

## Definition

Abe & Rajagopal's Tsallis conditional entropy between discrete random variables
``X`` and ``Y`` with finite ranges ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
H_q^{T_A}(X | Y) = \\dfrac{H_q^T(X, Y) - H_q^T(Y)}{1 + (1-q)H_q^T(Y)},
```

where ``H_q^T(\\cdot)`` and ``H_q^T(\\cdot, \\cdot)`` is the [`Tsallis`](@ref)
entropy and the joint Tsallis entropy.
"""
Base.@kwdef struct ConditionalEntropyTsallisAbe{B, Q} <: ConditionalEntropy
    base::B = 2
    q::Q = 1.5
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:ConditionalEntropyTsallisAbe}, inputs...)
    probs = probabilities(est.discretization, inputs...)
    return association(est.definition, probs)
end

function association(definition::ConditionalEntropyTsallisAbe, pxy::Probabilities{T, 2}) where {T}
    (; base, q) = definition

    if q == 1 # if shannon, normalize
        return association(ConditionalEntropyShannon(; base), pxy)
    end

    py = marginal(pxy, dims = 2)
    # Definition 7 in Abe & Rajagopal (2001)
    hjoint = 1 / (1 - q) * (sum(pxy .^ 2) - 1)

    # The marginal Tsallis entropy for the second variable
    hy = information(Tsallis(; q, base), py)

    # Equation 13 in Abe & Rajagopal (2001)
    ce = (hjoint - hy) / (1 + (1 - q)*hy)

    return ce
end