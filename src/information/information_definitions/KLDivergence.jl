export KLDivergence

"""
    KLDivergence <: BivariateInformationMeasure

The Kullback-Leibler (KL) divergence.

## Description

The KL-divergence between two probability distributions
``P_X = (p_x(\\omega_1), \\ldots, p_x(\\omega_n))`` and
``P_Y = (p_y(\\omega_1), \\ldots, p_y(\\omega_m))``, both defined over the same
[`OutcomeSpace`](@ref) ``\\Omega = \\{\\omega_1, \\ldots, \\omega_n \\}``, is defined as

```math
D_{KL}(P_Y(\\Omega) || P_Y(\\Omega)) =
\\sum_{\\omega \\in \\Omega} p_x(\\omega) \\log\\dfrac{p_x(\\omega)}{p_y(\\omega)}
```

## Implements

- [`information`](@ref). Used to compute the KL-divergence between two pre-computed
    probability distributions. If used with [`RelativeAmount`](@ref), the KL divergence may
    be undefined to due some outcomes having zero counts. Use some other
    [`ProbabilitiesEstimator`](@ref) like [`BayesianRegularization`](@ref) to ensure
    all estimated probabilities are nonzero.
"""
struct KLDivergence{B} <: BivariateInformationMeasure
    base::B
end
KLDivergence(; base = 2) = KLDivergence(base)

function information(measure::KLDivergence, px::Probabilities, py::Probabilities)
    size_match(measure, px, py)
    return sum(pxᵢ * log(measure.base, pxᵢ / pyᵢ) for (pxᵢ, pyᵢ) in zip(px, py))
end
