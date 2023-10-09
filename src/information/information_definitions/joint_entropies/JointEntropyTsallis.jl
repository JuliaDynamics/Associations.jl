using ComplexityMeasures: Tsallis

export JointEntropyTsallis

"""
    JointEntropyTsallis <: JointEntropy
    JointEntropyTsallis(; base = 2, q = 1.5)

The Tsallis joint entropy definition from [Furuichi2006](@citet). 

## Definition

Given two two discrete random variables ``X`` and ``Y`` with ranges ``\\mathcal{X}`` and
``\\mathcal{X}``, [Furuichi2006](@citet) defines the Tsallis joint entropy as

```math
H_q^T(X, Y) = -\\sum_{x\\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y)^q \\log_q p(x, y)
```

where ``log_q(x, q) = \\dfrac{x^{1-q} - 1}{1-q}`` is the q-logarithm.
"""
Base.@kwdef struct JointEntropyTsallis{B, Q} <: JointEntropy
    base::B = 2
    q::Q = 1.5
end

function information(definition::JointEntropyTsallis, pxy::Probabilities{T, 2}) where T
    (; base, q) = definition
    
    h = 0.0
    for p in pxy
        if p != 0.0 # logq will error if probability is zero.
            h += p^q * logq(p, q)
        end
    end
    h = -h
    return _convert_logunit(h, ℯ, base)
end