export CrossDifferentialEntropyEstimator
export entropy_cross

"""
    CrossEntropy <: InformationMeasure

The supertype of all cross-entropies.
"""
abstract type CrossEntropy <: InformationMeasure end

"""
The supertype of all cross-entropy definitions.
"""
abstract type CrossEntropyDefinition <: Definition end

"""
    CrossDifferentialEntropyEstimator <: InformationEstimator

The supertype of all cross-entropy estimators.

## Subtypes

- [`BulinskiDimitrov`](@ref).
"""
abstract type CrossDifferentialEntropyEstimator end

"""
    entropy_cross(measure::CrossEntropy, est::ProbabilitiesEstimator, x, y)

    entropy_cross(measure::CrossEntropy, est::CrossDifferentialEntropyEstimator, x, y)

Estimate `measure`, a cross-entropy of some kind, between `x` and `y` using the given
estimator.

- The first set of signatures is for discrete cross-entropies (meaning that they're
    computed directly from a set of probabilities estimated by a
    [`ProbabilitiesEstimator`](@ref)).
    For these methods to give meaningful results, you must ensure that the `est` gives
    the *the same* [`outcome_space`](@ref) for both `x` and `y`. For example, use a
    [`ValueHistogram`](@ref) with a [`FixedRectangularBinning`](@ref).
- The second set of signatures is for differential cross-entropies (meaning that they're
    computed from density estimates). For a full list of compatible definitions and
    estimators, see the online documentation.

Returns `div`, the divergence estimate, whose interpretation depends on the
combination of `definition` and `est`.

## Supported definitions

Cross-entropies appear in several forms in the literature. We currently support
the following measures:

- [`CrossEntropyShannon`](@ref). Discrete Shannon relative entropy.
"""
function entropy_cross end

entropy_cross(measure::CrossEntropy, est, x, y) = estimate(measure, est, x, y)
entropy_cross(def::CrossEntropyDefinition, measure::CrossEntropy, est, x, y) =
    estimate(def, measure, est, x, y)

include("measures/CrossEntropyShannon.jl")
include("measures/CrossEntropyRenyi.jl")

include("definitions/Thierrin.jl")

include("estimators/BulinskiDimitrov.jl")

# """
#     entropy_cross(e::EntropyDefinition, est::CrossDifferentialEntropyEstimator, x::AbstractDataset{D},
#         y::AbstractDataset{D}) where D

# Estimate the differential cross-entropy of type `e` between `D`-dimensional dataset
# `x` and `y`, using the provided [`CrossDifferentialEntropyEstimator`](@ref).

# ## Description

# The definition of cross entropy varies depending on which [`Entropy`](@ref) is it based on.
# Individual [`CrossDifferentialEntropyEstimator`](@ref)s specify the formulas they approximate.

# """
# function entropy_cross()

# end

# """
#     divergence(::Renyi, p::Probabilities, q::Probabilities)

# Estimate the (discrete) relative entropy, or KL divergence, between the pre-computed
# probability distributions `p` and `q`, where `p[i]` and `q[i]` is the probability of the
# `i`-th outcome in some [outcome_space](@ref) ``\\omega{X}``, defined as

# ```math
# D_{KL}(X || Y) = \\sum_{x \\in \\mathcal{X}} P(x) \\log{\\dfrac{P(x)}{Q(x)}}
# ```

# See also: [`probabilities`](@ref).
# """

# TODO: we can also define C(ℙ, ℚ) = D(ℙ || ℚ) + H(ℙ), so cross entropy
# can be estimated using any combination of relative entropy and entropy estimators,
# if the estimators are both defined for the same type of entropy (is this true? or
# only for Shannon entropy?)
# Is it wise to mix estimators in this way? Probably not always, but it is technically
# possible.
#include("estimators/estimators.jl")
