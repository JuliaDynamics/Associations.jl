export RelativeEntropyShannon
export ShannonDivergence

"""
    ShannonDivergence  <: DivergenceDefinition
    ShannonDivergence

An instruction to compute the discrete Shannon relative entropy (KL divergence) according to
the definition in ...
"""
struct ShannonDivergence <: DivergenceDefinition end

"""
    RelativeEntropyShannon <: Divergence
    RelativeEntropyShannon(e::Entropy = Shannon(), definition = ShannonDivergence())
    RelativeEntropyShannon(; base::Real)

`RelativeEntropyShannon` is a directive to compute the discrete Shannon relative entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

- [`ShannonDivergence`](@ref).

See also: [`divergence`](@ref).
"""
struct RelativeEntropyShannon{D <: DivergenceDefinition, E <: Entropy} <: Divergence
    e::E
    definition::D
    function RelativeEntropyShannon(; base = 2, q = 1.5,
            definition::D = ShannonDivergence()) where {D}
            e = Shannon(; base, q)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(measure::RelativeEntropyShannon{<:ShannonDivergence},
        est::ProbabilitiesEstimator, x, y)
    q, base = measure.e.q, measure.e.base
    p = probabilities(est, x)
    q = probabilities(est, y)
    re = 1 / (q - 1) * log(sum(pᵢ^q * qᵢ^(1 - q) for (pᵢ, qᵢ) in zip(p, q)))
    return re / log(base, ℯ)
end
