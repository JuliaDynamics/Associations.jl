export RenyiDivergence

"""
    RenyiDivergence <: BivariateInformationMeasure
    RenyiDivergence(q; base = 2)

The Rényi divergence of positive order `q`

## Description

The Rényi divergence between two probability distributions
``P_X = (p_x(\\omega_1), \\ldots, p_x(\\omega_n))`` and
``P_Y = (p_y(\\omega_1), \\ldots, p_y(\\omega_m))``, both defined over the same
[`OutcomeSpace`](@ref) ``\\Omega = \\{\\omega_1, \\ldots, \\omega_n \\}``, is defined as
[vanErven2014](@cited)

```math
D_{q}(P_Y(\\Omega) || P_Y(\\Omega)) =
\\dfrac{1}{q - 1} \\log \\sum_{\\omega \\in \\Omega}p_x(\\omega)^{q}p_y(\\omega)^{1-\\alpha}
```

## Implements

- [`information`](@ref). Used to compute the Rényi divergence between two pre-computed
    probability distributions. If used with [`RelativeAmount`](@ref), the KL divergence may
    be undefined to due some outcomes having zero counts. Use some other
    [`ProbabilitiesEstimator`](@ref) like [`BayesianRegularization`](@ref) to ensure
    all estimated probabilities are nonzero.
"""
struct RenyiDivergence{Q, B} <: BivariateInformationMeasure
    q::Q
    base::B
    function RenyiDivergence(q::Q, base::B) where {Q, B}
        q > 0 || throw(ArgumentError("`q` must be positive. Got $q"))
        new{Q, B}(q, base)
    end
end
RenyiDivergence(;q = 0.5, base = 2) = RenyiDivergence(q, base)

function information(measure::RenyiDivergence, px::Probabilities, py::Probabilities)
    base, q = measure.base, measure.q

    if q == Inf
        return maximum(pxᵢ / pyᵢ for (pxᵢ, pyᵢ) in zip(px, py))
    end
    s = 0.0
    for (pxᵢ, pyᵢ) in zip(px, py)
        if pxᵢ != 0.0 && pyᵢ != 0.0
            s += pxᵢ^q * pyᵢ^(1 - q)
        end
    end
    return 1 / (q - 1) * log(base, s)
end
