export RelativeEntropyEstimator
export entropy_relative

abstract type RelativeEntropyEstimator end

#include("tsallis/entropy_relative_tsallis.jl")

"""
    entropy_relative([e::Entropy], est::RelativeEntropyEstimator,
        x::Vector_or_Dataset, y::Vector_or_Dataset)

Compute the (differential) relative entropy of type `e` between `x` and `y` using the
provided [`RelativeEntropyEstimator`](@ref).

The first argument, the entropy type, is optional; it defaults to [`Shannon`](@ref).
"""
function entropy_relative end

include("estimators/Wang2009.jl")
include("estimators/BulinskiDimitrov.jl")
