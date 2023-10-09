using ComplexityMeasures: Shannon

export JointEntropyShannon

"""
    JointEntropyShannon <: JointEntropy
    JointEntropyShannon(; base = 2)

The Shannon joint entropy measure [CoverThomas2006](@cite).

## Definition

Given two two discrete random variables ``X`` and ``Y`` with ranges ``\\mathcal{X}`` and
``\\mathcal{X}``, [CoverThomas2006](@citet) defines the Shannon joint entropy as

```math
H^S(X, Y) = -\\sum_{x\\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y) \\log p(x, y)
```
"""
Base.@kwdef struct JointEntropyShannon{B} <: JointEntropy
    base::B = 2
end

function information(definition::JointEntropyShannon, pxy::Probabilities{T, 2}) where T
    (; base) = definition
    
    h = 0.0
    for p in pxy
        h += p * log(p)
    end
    h = -h
    return _convert_logunit(h, ℯ, base)
end