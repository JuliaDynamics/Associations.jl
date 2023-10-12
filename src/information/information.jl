import ComplexityMeasures: information 

include("core.jl")

# These files extend the single-variable information API in ComplexityMeasures.jl.
include("information_functions.jl")
include("information_estimators/information_estimators.jl")
include("information_definitions/information_definitions.jl")

# Specific estimators must be included after definitions.

include("information_estimators/mutual_info_estimators/mutual_info_estimators.jl")
include("information_estimators/conditional_mutual_info_estimators/conditional_mutual_info_estimators.jl")

include("convenience.jl")