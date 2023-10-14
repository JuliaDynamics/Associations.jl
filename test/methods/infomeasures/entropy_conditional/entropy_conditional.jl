@test ConditionalEntropyShannon() isa ConditionalEntropyShannon
@test ConditionalEntropyTsallisAbe() isa ConditionalEntropyTsallisAbe
@test ConditionalEntropyTsallisFuruichi() isa ConditionalEntropyTsallisFuruichi

include("contingency_matrix.jl")
include("discrete.jl")
include("continuous.jl")
