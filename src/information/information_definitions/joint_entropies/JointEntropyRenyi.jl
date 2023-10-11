using ComplexityMeasures: Renyi

export JointEntropyRenyi

"""
    JointEntropyRenyi <: JointEntropy
    JointEntropyRenyi(; base = 2, q = 1.5)

The Rényi joint entropy measure [Golshani2009](@cite).

## Definition

Given two two discrete random variables ``X`` and ``Y`` with ranges ``\\mathcal{X}`` and
``\\mathcal{X}``, [Golshani2009](@citet) defines the Rényi joint entropy as

```math
H_q^R(X, Y) = \\dfrac{1}{1-\\alpha} \\log \\sum_{i = 1}^N p_i^q,
```

where ``q > 0`` and ``q != 1``.
"""
Base.@kwdef struct JointEntropyRenyi{B, Q} <: JointEntropy
    base::B = 2
    q::Q = 1.5
end

function information(definition::JointEntropyRenyi, pxy::Probabilities{T, 2}) where T
    (; base, q) = definition
    
    h = 0.0
    for p in pxy
        if p != 0
            h += p^q
        end
    end
    h = 1 / (1 - q) * log(h) 
    return _convert_logunit(h, ℯ, base)
end