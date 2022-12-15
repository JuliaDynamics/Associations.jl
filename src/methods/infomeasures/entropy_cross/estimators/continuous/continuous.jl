# """
#     ContinuousCrossEntropy <: CrossEntropyEstimator
#     ContinuousCrossEntropy(est::ProbabilitiesEstimator, definition::CrossEntropyDefinition)

# `ContinuousCrossEntropy` is a generic plug-in estimator for continuous generalized cross entropy.

# ## Compatible definitions

# Multiple generalizations of cross entropy exists, and `definition` specifies which
# generalization is to be computed. Currently, we support:

# ## Description

# The `ContinuousCrossEntropy estimator first uses `est` to compute probabilities, then plugs
# these probabilities into the cross entropy formula specified by `definition`.
# """
# struct ContinuousCrossEntropy{P, D <: CrossEntropyDefinition} <: CrossEntropyEstimator
#     est::P
#     definition::D
# end

include("BulinskiDimitrov.jl")
include("RenyiCrossEntropy.jl")
