export VariationDistance

"""
    VariationDistance <: BivariateInformationMeasure

The Variation distance.

## Description

The Variation distance between two probability distributions
``P_X = (p_x(\\omega_1), \\ldots, p_x(\\omega_n))`` and
``P_Y = (p_y(\\omega_1), \\ldots, p_y(\\omega_m))``, both defined over the same
[`OutcomeSpace`](@ref) ``\\Omega = \\{\\omega_1, \\ldots, \\omega_n \\}``, is
[defined](https://en.wikipedia.org/wiki/Variation_distance) as

```math
D_{V}(P_Y(\\Omega) || P_Y(\\Omega)) =
\\dfrac{1}{2} \\sum_{\\omega \\in \\Omega} | p_x(\\omega) - p_y(\\omega) |
```

## Implements

- [`information`](@ref). Used to compute the variation distance between two pre-computed
    probability distributions.
"""
struct VariationDistance <: BivariateInformationMeasure end

function information(measure::VariationDistance, px::Probabilities, py::Probabilities)
    return 1/2 * sum(abs(pxᵢ - pyᵢ) for (pxᵢ, pyᵢ) in zip(px, py))
end
