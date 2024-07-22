using ComplexityMeasures: Tsallis

export JointEntropyTsallis

"""
    JointEntropyTsallis <: JointEntropy
    JointEntropyTsallis(; base = 2, q = 1.5)

The Tsallis joint entropy definition from [Furuichi2006](@citet). 

## Usage 

- Use with [`association`](@ref) to compute the Furuichi-Tsallis joint entropy between 
    two variables.

## Compatible estimators

- [`JointProbabilities`](@ref)

## Definition

Given two two discrete random variables ``X`` and ``Y`` with ranges ``\\mathcal{X}`` and
``\\mathcal{X}``, [Furuichi2006](@citet) defines the Tsallis joint entropy as

```math
H_q^T(X, Y) = -\\sum_{x\\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y)^q \\log_q p(x, y),
```

where ``log_q(x, q) = \\dfrac{x^{1-q} - 1}{1-q}`` is the q-logarithm, and 
we define ``log_q(x, q) := 0`` if ``q = 0``.

## Estimation

- [Example 1](@ref example_JointEntropyTsallis_OrdinalPatterns): 
    [`JointProbabilities`](@ref) with [`OrdinalPatterns`](@ref) outcome space
"""
Base.@kwdef struct JointEntropyTsallis{B, Q} <: JointEntropy
    base::B = 2
    q::Q = 1.5
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:JointEntropyTsallis}, x, y)
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

function association(definition::JointEntropyTsallis, pxy::Probabilities{T, 2}) where T
    (; base, q) = definition
    
    h = 0.0
    for p in pxy
        if p != 0.0 # Define logq(0) = 0
            h += p^q * logq(p, q)
        end
    end
    h = -h
    return _convert_logunit(h, ℯ, base)
end