"""
The supertype for divergence estimators.

Currently implemented subtypes are:

- [`PoczosSchneiderRE`](@ref)
- [`Wang`](@ref)
- [`WangTransformed`](@ref)
"""
abstract type DivergenceEstimator <: InformationMeasureEstimator end

include("renyi/RelativeEntropyRenyi.jl")
include("renyi/RelativeEntropyRenyiDifferential.jl")
include("tsallis/RelativeEntropyTsallis.jl")
include("tsallis/RelativeEntropyTsallisDifferential.jl")
include("shannon/RelativeEntropyShannon.jl")
include("shannon/RelativeEntropyShannonDifferential.jl")

include("dedicated_estimators/estimators.jl")
