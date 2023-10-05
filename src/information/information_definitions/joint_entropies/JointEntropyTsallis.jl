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
struct JointEntropyTsallis{E<:Tsallis} <: JointEntropy
    e::E
    function JointEntropyTsallis(; base = 2, q = 1.5)
        e = Tsallis(; base, q)
        new{typeof(e)}(e)
    end
end

function information(definition::JointEntropyTsallis, pxy::Probabilities{T, 2}) where T
    base = definition.e.base
    q = definition.e.q
    h = 0.0
    for p in pxy
        h += p^q * logq(p, q)
    end
    h = -h
    return _convert_logunit(h, ℯ, base)
end