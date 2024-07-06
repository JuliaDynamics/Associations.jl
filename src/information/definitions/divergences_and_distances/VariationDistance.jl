export VariationDistance

"""
    VariationDistance <: DivergenceOrDistance

The variation distance.

## Usage 

- Use with [`association`](@ref) to compute the compute the variation distance between two 
    pre-computed probability distributions, or from raw data using one of the estimators
    listed below.

## Compatible estimators

- [`JointDistanceDistribution`](@ref)

## Description

The variation distance between two probability distributions
``P_X = (p_x(\\omega_1), \\ldots, p_x(\\omega_n))`` and
``P_Y = (p_y(\\omega_1), \\ldots, p_y(\\omega_m))``, both defined over the same
[`OutcomeSpace`](@ref) ``\\Omega = \\{\\omega_1, \\ldots, \\omega_n \\}``, is
[defined](https://en.wikipedia.org/wiki/Variation_distance) as

```math
D_{V}(P_Y(\\Omega) || P_Y(\\Omega)) =
\\dfrac{1}{2} \\sum_{\\omega \\in \\Omega} | p_x(\\omega) - p_y(\\omega) |
```

## Examples

```julia
using CausalityTools
# From raw data
using Random; rng = Xoshiro(1234)
n = 100000
x, y = rand(rng, n), rand(rng, n)
est = JointProbabilities(VariationDistance(), CodifyVariables(OrdinalPatterns(m=3)))
div_hd = association(est, x, y) # pretty close to zero

# From pre-computed PMFs
p1 = Probabilities([0.1, 0.5, 0.2, 0.2])
p2 = Probabilities([0.3, 0.3, 0.2, 0.2])
association(VariationDistance(), p1, p2)
```
"""
struct VariationDistance <: DivergenceOrDistance end

function information(measure::VariationDistance, px::Probabilities, py::Probabilities)
    return 1/2 * sum(abs(pxᵢ - pyᵢ) for (pxᵢ, pyᵢ) in zip(px, py))
end
