abstract type DivergenceOrDistance <: BivariateInformationMeasure end

include("KLDivergence.jl")
include("RenyiDivergence.jl")
include("HellingerDistance.jl")
include("VariationDistance.jl")

