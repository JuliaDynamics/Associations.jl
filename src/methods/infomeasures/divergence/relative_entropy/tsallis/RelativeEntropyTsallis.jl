export RelativeEntropyTsallis
export TsallisDivergence

"""
    TsallisDivergence <: DivergenceDefinition
    TsallisDivergence()

`TsallisDivergenceTsallis1998` gives the original definition of Tsallis relative entropy
from Tsallis (1998)[^Tsallis1998].

## Description

For a discrete sample space ``\\Omega`` and probability mass functions
``a(x) : \\Omega \\to [0, 1]`` and ``b(x) : \\Omega \\to [0, 1]``,
the Tsallis relative entropy (divergence) is given by

```math
D_q^T(A || B) = - \\sum_{i = 1}^n a_i^q \\log_q \\dfrac{b_i}{a_i}.
```

For `q = 1`, the Shannon relative entropy (KL divergence) is obtained, and ``log_q`` is
the q-logarithm defined in Tsallis (1998).

[^Tsallis1998]:
    Tsallis, C. (1998). Generalized entropy-based criterion for consistent testing.
    Physical Review E, 58(2), 1442.
"""
struct TsallisDivergence <: DivergenceDefinition end

"""
    RelativeEntropyTsallis <: Divergence
    RelativeEntropyTsallis(e::Entropy = Shannon(), definition = TsallisDivergence())
    RelativeEntropyTsallis(; base::Real)

`RelativeEntropyTsallis` is a directive to compute the discrete Tsallis relative entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

- [`TsallisDivergence`](@ref).

See also: [`divergence`](@ref).
"""
struct RelativeEntropyTsallis{D <: DivergenceDefinition, E <: Entropy} <: Divergence
    e::E
    definition::D
    function RelativeEntropyTsallis(; base = 2, q = 1.5,
            definition::D = TsallisDivergence()) where {D}
            e = Tsallis(; base, q)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(measure::RelativeEntropyTsallis{<:TsallisDivergence},
        est::ProbabilitiesEstimator, x, y)
    q, base = measure.e.q, measure.e.base
    P = probabilities(est, x)
    Q = probabilities(est, y)
    return -sum(pᵢ * qlog(qᵢ / pᵢ, q; base) for (pᵢ, qᵢ) in zip(P, Q))
end

function qlog(x, q; base = 2)
    if x <= 0
        return 0
    else
        if q == 1 # base only affects result if q == 1 (i.e. we get Shannon divergence)
            return Entropies.log_with_base(base)(x)
        else
            return (x^(1 - q) - 1) / (1 - q)
        end
    end
end
