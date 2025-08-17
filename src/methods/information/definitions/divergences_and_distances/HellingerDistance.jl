export HellingerDistance

"""
    HellingerDistance <: DivergenceOrDistance

The Hellinger distance.

## Usage 

- Use with [`association`](@ref) to compute the compute the Hellinger distance between two pre-computed
    probability distributions, or from raw data using one of the estimators listed below.

## Compatible estimators

- [`JointProbabilities`](@ref)

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

## Estimation

- [Example 1](@ref example_HellingerDistance_precomputed_probabilities): From precomputed probabilities
- [Example 2](@ref example_HellingerDistance_JointProbabilities_OrdinalPatterns): 
    [`JointProbabilities`](@ref) with [`OrdinalPatterns`](@extref ComplexityMeasures) outcome space
"""
struct HellingerDistance <: DivergenceOrDistance end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:HellingerDistance}, x, y)
    # Dispatch to generic method in `divergences_and_distances.jl` with 2D `Probabilities`
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

function association(measure::HellingerDistance, px::Probabilities, py::Probabilities)
    return 1 / sqrt(2) * sum((sqrt(pxᵢ) - sqrt(pyᵢ))^2 for (pxᵢ, pyᵢ) in zip(px, py))
end
