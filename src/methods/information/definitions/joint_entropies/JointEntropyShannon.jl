using ComplexityMeasures: Shannon

export JointEntropyShannon

"""
    JointEntropyShannon <: JointEntropy
    JointEntropyShannon(; base = 2)

The Shannon joint entropy measure [CoverThomas1999](@cite).

## Usage 

- Use with [`association`](@ref) to compute the Shannon joint entropy between 
    two variables.

## Compatible estimators

- [`JointProbabilities`](@ref)

## Definition

Given two two discrete random variables ``X`` and ``Y`` with ranges ``\\mathcal{X}`` and
``\\mathcal{X}``, [CoverThomas1999](@citet) defines the Shannon joint entropy as

```math
H^S(X, Y) = -\\sum_{x\\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y) \\log p(x, y),
```

where we define ``log(p(x, y)) := 0`` if ``p(x, y) = 0``.

## Estimation

- [Example 1](@ref example_JointEntropyShannon_Dispersion): 
    [`JointProbabilities`](@ref) with [`Dispersion`](@ref) outcome space
"""
Base.@kwdef struct JointEntropyShannon{B} <: JointEntropy
    base::B = 2
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:JointEntropyShannon}, x, y)
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

function association(definition::JointEntropyShannon, pxy::Probabilities{T, 2}) where T
    (; base) = definition
    
    h = 0.0
    for p in pxy
        if p != 0 # Define log(0) = 0
            h += p * log(p)
        end
    end
    h = -h
    return _convert_logunit(h, ℯ, base)
end