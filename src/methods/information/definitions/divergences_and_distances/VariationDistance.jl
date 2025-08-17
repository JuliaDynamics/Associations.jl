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

- [Example 1](@ref example_VariationDistance_precomputed_probabilities): From precomputed probabilities
- [Example 2](@ref example_VariationDistance_JointProbabilities_OrdinalPatterns): 
    [`JointProbabilities`](@ref) with [`OrdinalPatterns`](@extref ComplexityMeasures) outcome space
"""
struct VariationDistance <: DivergenceOrDistance end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:VariationDistance}, x, y)
    # Dispatch to generic method in `divergences_and_distances.jl` with 2D `Probabilities`
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

function association(measure::VariationDistance, px::Probabilities, py::Probabilities)
    return 1 / 2 * sum(abs(pxᵢ - pyᵢ) for (pxᵢ, pyᵢ) in zip(px, py))
end
