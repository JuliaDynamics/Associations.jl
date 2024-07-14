include("counts_and_probs/counts_and_probs.jl")
include("core.jl")

# These files extend the single-variable information API in ComplexityMeasures.jl.
include("information_functions.jl")
include("estimators/information_estimators.jl")
include("definitions/information_definitions.jl")

# Specific estimators must be included after definitions.
include("estimators/mutual_info_estimators/mutual_info_estimators.jl")
include("estimators/conditional_mutual_info_estimators/conditional_mutual_info_estimators.jl")
include("estimators/transfer_entropy_estimators/transfer_entropy_estimators.jl")

include("convenience.jl")