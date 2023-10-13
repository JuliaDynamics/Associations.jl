include("discrete/joint.jl")
include("entropy_decomposition.jl")

# Generically test the estimators. Any future analytical tests should go in their own 
# estimator-specific files.
include("cond_mutual_information_estimators/cond_mutual_information_estimators.jl")
include("mutual_information_estimators/mutual_information_estimators.jl")