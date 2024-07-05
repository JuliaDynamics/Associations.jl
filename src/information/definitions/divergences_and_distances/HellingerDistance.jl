export HellingerDistance

"""
    HellingerDistance <: DivergenceOrDistance

The Hellinger distance.

## Usage 

- [`information`](@ref). Used to compute the Hellinger distance between two pre-computed
    probability distributions.

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

## Examples

```julia
# There should be zero information gain from `x` over `y` for independent random variables.
using Random; rng = Xoshiro(1234)
n = 100000
x, y = rand(rng, n), rand(rng, n)
o = OrdinalPatterns(m = 3)
div_hd = information(HellingerDistance(), o, x, y) # pretty close to zero
```

"""
struct HellingerDistance <: DivergenceOrDistance end

function information(measure::HellingerDistance, px::Probabilities, py::Probabilities)
    return 1/sqrt(2) * sum((sqrt(pxᵢ) - sqrt(pyᵢ))^2 for (pxᵢ, pyᵢ) in zip(px, py))
end
