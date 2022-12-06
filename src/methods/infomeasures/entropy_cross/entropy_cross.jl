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
    entropy_cross([e::Entropy,] est::CrossEntropyEstimator,
        x::AbstractDataset{D}, y::AbstractDataset{D})

Estimate the differential cross-entropy of type `e` between `x` and `y`, using
the provided [`CrossEntropyEstimator`](@ref).

The first argument, the entropy type, is optional; it defaults to `Shannon(; base = 2)`.

## Description

## [`Shannon`](@ref) cross entropy

Following the notation of Bulinski & Dimitrov (2021)[^Bulinski2021],
let ``\\mathbb{P}`` and ``\\mathbb{Q}`` be continuous probability measures
with densities ``p(x)`` and ``q(x)``, ``x \\in \\mathcal{R}^D``,
with respect to the Lebesque measure ``\\mu``. Then, writing ``dx := \\mu(dx)``,
Shannon cross entropy is defined as

```math
C(\\mathbb{P}, \\mathbb{Q}) = - \\int_{\\mathbb{R}^d} p(x) \\log{(q(x))} dx.
```

## Other cross entropies

Besides the Shannon cross entropy, one can also define cross entropies based on
other generalized entropies. However, there is no consensus in the literature on
what the definition of such cross  entropies are. Therefore, individual
[`CrossEntropyEstimator`](@ref) specify precisely which quantities they approximate.

See also: [`Entropy`](@ref).
"""
function entropy_cross(e::Entropy, est::CrossEntropyEstimator,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D
end

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
