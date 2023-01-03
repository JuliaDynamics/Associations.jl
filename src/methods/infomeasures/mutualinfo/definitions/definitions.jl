
# There are multiple ways of formulating the Shannon mutual information and estimating it.
# Each struct below defines a unique way of doing so, and implements relevant methods
# for `estimate`, if not using generic dispatch as defined in `common_dispatch.jl`.
include("MIDefinitionShannonDoubleSum.jl")
include("MIDefinitionShannonH3.jl")
# TODO: new estimation using KL divergence when implemented

# There are multiple ways of formulating the Tsallis mutual information and estimating it.
# Each file below defines a unique way of doing so, and implements relevant methods
# for `estimate`, if not using generic dispatch as defined in `common_dispatch.jl`.
include("MIDefinitionTsallisH3Furuichi.jl")
include("MIDefinitionTsallisH3Martin.jl")
# TODO: new estimation using KL divergence when implemented
# A different version of Tsallis MI is given in: https://www.mdpi.com/1099-4300/17/8/5382

include("MIDefinitionRenyiSarbu.jl")
