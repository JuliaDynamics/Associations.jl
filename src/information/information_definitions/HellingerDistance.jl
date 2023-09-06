export HellingerDistance

"""
    HellingerDistance <: BivariateInformationMeasure

The Hellinger distance.

## Description

The Hellinger distance between two probability distributions
``P_X = (p_x(\\omega_1), \\ldots, p_x(\\omega_n))`` and
``P_Y = (p_y(\\omega_1), \\ldots, p_y(\\omega_m))``, both defined over the same
[`OutcomeSpace`](@ref) ``\\Omega = \\{\\omega_1, \\ldots, \\omega_n \\}``, is
[defined](https://en.wikipedia.org/wiki/Hellinger_distance) as

```math
D_{H}(P_Y(\\Omega) || P_Y(\\Omega)) =
\\dfrac{1}{\\sqrt{2}} \\sum_{\\omega \\in \\Omega} (\\sqrt{p_x(\\omega)} - \\sqrt{p_y(\\omega)})^2
```

## Implements

- [`information`](@ref). If used with [`RelativeAmount`](@ref), the KL divergence may
    be undefined to due some outcomes having zero counts. Use some other
    [`ProbabilitiesEstimator`](@ref) like [`BayesianRegularization`](@ref) to ensure
    all estimated probabilities are nonzero.
## Implements

- [`information`](@ref). If used with [`RelativeAmount`](@ref), the KL divergence may
    be undefined to due some outcomes having zero counts. Use some other
    [`ProbabilitiesEstimator`](@ref) like [`BayesianRegularization`](@ref) to ensure
    all estimated probabilities are nonzero.
"""
struct HellingerDistance <: BivariateInformationMeasure end

function information(measure::HellingerDistance, px::Probabilities, py::Probabilities)
    return 1/sqrt(2) * sum((sqrt(pxᵢ) - sqrt(pyᵢ))^2 for (pxᵢ, pyᵢ) in zip(px, py))
end
