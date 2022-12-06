export entropy_conditional
export ConditionalEntropyEstimator

# There are at least three ways of defining conditional renyi entropy.
# https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6191351
# The Jizba estimator gives #2.
"""
The supertype for all conditional entropy estimators.
"""
abstract type ConditionalEntropyEstimator end

"""
    entropy_conditional(e::Entropy,
        est::ConditionalEntropyEstimator{ProbabilitiesEstimator}, x, y)

Compute the conditional entropy of type `e` of `x` given `y` using the given estimator.

The input data `x` and `y` can be either univariate vectors or multivariate datasets,
depending on `est`.
"""
function entropy_conditional end

include("estimators/estimators.jl")
