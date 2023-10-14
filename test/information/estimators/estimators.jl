include("discrete/joint.jl")

# Generically test the estimators. Any future analytical tests should go in their own 
# estimator-specific files.
include("transfer_entropy_estimators/transfer_entropy_estimators.jl")