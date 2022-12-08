export CrossEntropyEstimator
export entropy_cross

"""
    CrossEntropy <: InformationMeasure

The cross-entropy. Used with [`estimate`](@ref) to compute cross-entropy.
"""
struct CrossEntropy <: InformationMeasure end

"""
    CrossEntropyEstimator <: InformationEstimator

Subtypes of `CrossEntropyEstimator` estimate the differential cross-entropy
(see [`entropy_cross`](@ref)).

Currently implemented subtypes are:

- [`BulinskiDimitrov`](@ref).
"""
abstract type CrossEntropyEstimator end

"""
    entropy_cross(e::Entropy, est::CrossEntropyEstimator,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D

Estimate the differential cross-entropy of type `e` between `D`-dimensional dataset
 `x` and `y`, using the provided [`CrossEntropyEstimator`](@ref).

## Description

The definition of cross entropy varies depending on which [`Entropy`](@ref) is it based on.
Individual [`CrossEntropyEstimator`](@ref)s specify the formulas they approximate.
"""
function entropy_cross end

# """
#     entropy_relative(::Renyi, p::Probabilities, q::Probabilities)

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
include("estimators/BulinskiDimitrov.jl")
