"""
    DivergenceOrDistance <: BivariateInformationMeasure

The supertype for bivariate information measures aiming to quantify some sort of
divergence, distance or closeness between two probability distributions.

Some of these measures are proper metrics, while others are not, but they have in
common that they aim to quantify how "far from each other" two probabilities distributions
are.

## Concrete implementations

- [`HellingerDistance`](@ref)
- [`KLDivergence`](@ref)
- [`RenyiDivergence`](@ref)
- [`VariationDistance`](@ref)
"""
abstract type DivergenceOrDistance <: BivariateInformationMeasure end

# If a joint probability is given, get the marginals
function information(measure::DivergenceOrDistance, p::Probabilities{T, 2}) where T
    px = marginal(p, dims = 1)
    py = marginal(p, dims = 2)
    return information(measure, px, py)
end

include("KLDivergence.jl")
include("RenyiDivergence.jl")
include("HellingerDistance.jl")
include("VariationDistance.jl")

