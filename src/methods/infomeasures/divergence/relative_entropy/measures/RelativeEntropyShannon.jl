export RelativeEntropyShannon

"""
    RelativeEntropyShannon <: Divergence
    RelativeEntropyShannon(; base::Real)

`RelativeEntropyShannon` is a directive to compute the discrete Shannon relative entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

- [`ShannonDivergence`](@ref).

See also: [`divergence`](@ref).
"""
struct RelativeEntropyShannon{E <: EntropyDefinition} <: Divergence
    e::E
    function RelativeEntropyShannon(; base = 2)
            e = Shannon(; base)
        new{typeof(e)}(e)
    end
end

estimate(measure::RelativeEntropyShannon, est, x, y) =
    estimate(ShannonDivergence(), measure::RelativeEntropyShannon, est, x, y)

# function estimate(def::ShannonDivergence, measure::RelativeEntropyShannon,
#         est::ProbabilitiesEstimator, x, y)
#     base = measure.e.base
#     P = probabilities(est, x)
#     Q = probabilities(est, y)
#     re = sum(pᵢ * log0(base, pᵢ/qᵢ) for (pᵢ, qᵢ) in zip(P, Q))
#     return re / log(base, ℯ)
# end
