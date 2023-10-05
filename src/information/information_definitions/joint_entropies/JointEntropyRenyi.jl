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
H_q^R(X, Y) = -\\dfrac{1}{1-\\alpha} \\log \\sum_{i = 1}^N p_i^q.
```
"""
struct JointEntropyRenyi{E<:Renyi} <: JointEntropy
    e::E
    function JointEntropyRenyi(; base = 2, q = 1.5)
        e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
end

function information(definition::JointEntropyRenyi, pxy::Probabilities{T, 2}) where T
    base = definition.e.base
    q = definition.e.q
    h = 0.0
    for p in pxy
        h += p^q * logq(p, q)
    end
    h = -h
    return _convert_logunit(h, ℯ, base)
end