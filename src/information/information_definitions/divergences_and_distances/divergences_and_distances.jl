"""
    DivergenceOrDistance <: BivariateInformationMeasure

An abstract type for divergence or distance measures.

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

