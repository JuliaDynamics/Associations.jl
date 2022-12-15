"""
The supertype for divergence estimators.

Currently implemented subtypes are:

- [`PoczosSchneiderRE`](@ref)
- [`Wang`](@ref)
- [`WangTransformed`](@ref)
"""
abstract type DivergenceEstimator <: InformationMeasureEstimator end

include("definitions/ShannonDivergence.jl")
include("definitions/RenyiDivergence.jl")
include("definitions/TsallisDivergence.jl")
include("definitions/TsallisDivergenceDifferential.jl")

include("measures/RelativeEntropyRenyi.jl")
include("measures/RelativeEntropyTsallis.jl")
include("measures/RelativeEntropyShannon.jl")

include("estimators/estimators.jl")
