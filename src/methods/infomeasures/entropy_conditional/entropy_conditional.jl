export ConditionalEntropy
export ConditionalEntropyEstimator
export ConditionalEntropyDefinition
export entropy_conditional

# There are at least three ways of defining conditional renyi entropy.
# https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6191351
# The Jizba estimator gives #2.
""" The supertype of all conditional entropy measures """
abstract type ConditionalEntropy <: InformationMeasure end

""" The supertype of all conditional entropy definitions """
abstract type ConditionalEntropyDefinition <: Definition end

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

"""
    entropy_conditional(definition::ConditionalEntropyShannon, est::ProbabilitiesEstimator, x, y) → ce::Real
    entropy_conditional(definition::ConditionalEntropyTsallis, est::ProbabilitiesEstimator, x, y) → ce::Real

    entropy_conditional(definition::ConditionalEntropyShannonDifferential, est::EntropyEstimator, x, y) → ce::Real

Estimate the generalized conditional entropy (CE) (using the formula and logarithm base
specified by `definition`) between `x` and `y`, using the provided estimator.

The first set of signatures is for the discrete CE. The second set of signatures is for
differential/continuous CE For a full list of compatible definitions and estimators,
see the online documentation.

Returns `ce`, the conditional entropy estimate, whose interpretation depends on the
combination of `definition` and `est`.

## Supported definitions

Multiple generalized conditional entropies are found in the literature.
We currently support the following measures:

- [`ConditionalEntropyShannon`](@ref). Discrete Shannon conditional entropy.
- [`ConditionalEntropyTsallis`](@ref). Discrete Tsallis conditional entropy.
- [`ConditionalEntropyShannonDifferential`](@ref). Differential Shannon conditional entropy.
"""
entropy_conditional(measure, est, x, y) = estimate(measure, est, x, y)

include("shannon/ConditionalEntropyShannon.jl")
include("shannon/ConditionalEntropyShannonDifferential.jl")
include("tsallis/ConditionalEntropyTsallis.jl")
include("renyi/ConditionalEntropyRenyi.jl")
