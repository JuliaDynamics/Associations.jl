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

Compute the conditional entropy of type `e` of `x` given `y`, which we denote
``h^e(X|Y)`` in the continuous case and ``H^e(X|Y)`` in the discrete case, using the given
estimator.

The input data `x` and `y` can be either univariate vectors or multivariate datasets,
depending on `est`.

## Description

The definition of the conditional entropy varies depending on which [`Entropy`](@ref) is
it based on. Individual [`ConditionalEntropyEstimator`](@ref)s specify the formulas they
approximate.
"""
function entropy_conditional end

include("estimators/estimators.jl")
