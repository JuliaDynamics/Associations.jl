export RelativeEntropyTsallis

"""
    RelativeEntropyTsallis <: Divergence
    RelativeEntropyTsallis(; base = 2)

`RelativeEntropyTsallis` is a directive to compute the discrete Tsallis relative entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

- [`TsallisDivergence`](@ref).

See also: [`divergence`](@ref).
"""
struct RelativeEntropyTsallis{E <: EntropyDefinition} <: Divergence
    e::E
    function RelativeEntropyTsallis(; base = 2, q = 1.5)
            e = Tsallis(; base, q)
        new{typeof(e)}(e)
    end
end

estimate(measure::RelativeEntropyTsallis, est, x, y) =
    estimate(TsallisDivergence(), measure::RelativeEntropyTsallis, est, x, y)


function estimate(def::TsallisDivergence, measure::RelativeEntropyTsallis,
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
