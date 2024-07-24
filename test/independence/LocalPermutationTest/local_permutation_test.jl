include("api.jl")

# Measure-specific implementations. One file per method that is listed in the docstring 
# of `LocalPermutationTest`.
include("conditional_mutual_information.jl")
include("part_mutual_information.jl")
#include("transferentropy.jl")
include("partial_correlation.jl")
include("distance_correlation.jl")