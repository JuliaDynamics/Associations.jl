export RelativeEntropyEstimator
export entropy_relative

"""
The supertype for differential relative entropy (see [`entropy_relative`](@ref)) (also
called divergence) estimators.

Currently implemented subtypes are:

- [`PoczosSchneiderRE`](@ref)
- [`Wang`](@ref)
- [`WangTransformed`](@ref)
"""
abstract type RelativeEntropyEstimator end

"""
    entropy_relative(e::Entropy, est::RelativeEntropyEstimator,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D

Estimate the (differential) relative entropy of type `e` between the `D`-dimensional
datasets `x` and `y` using the provided [`RelativeEntropyEstimator`](@ref).

## Description

The definition of cross entropy varies depending on which [`Entropy`](@ref) is it based on.
Individual [`CrossEntropyEstimator`](@ref)s specify the formulas they approximate.

[^Bulinski2021]:
    Bulinski, A., & Dimitrov, D. (2021). Statistical estimation of the Kullback-Leibler
    divergence. Mathematics, 9(5), 544.
"""
function entropy_relative end
# entropy_relative(est::RelativeEntropyEstimator, args...; base = 2, kwargs...) =
#     entropy_relative(Shannon(; base), est, args...; kwargs...)


entropy_relative(e::Entropy, est::RelativeEntropyEstimator, x, y) =
    entropy_relative(e, est, Dataset(x), Dataset(y))
# """
#     entropy_relative(::Renyi, p::Probabilities, q::Probabilities)

# Estimate the (discrete) relative entropy, or KL divergence, between the pre-computed
# probability distributions `p` and `q`, where `p[i]` and `q[i]` is the probability of the
# `i`-th outcome in some [outcome_space](@ref) ``\\omega{X}``, defined as

# ```math
# D_{KL}(X || Y) = \\sum_{x \\in \\mathcal{X}} P(x) \\log{\\dfrac{P(x)}{Q(x)}}
# ```

# For this definition to be meaningful `p` and `q` must have *the same* outcome space.

# See also: [`probabilities`](@ref).
# """
function entropy_relative end

using Distributions: MvNormal
using LinearAlgebra: tr, inv
function entropy_relative(e::Renyi, N1::MvNormal, N2::MvNormal)
    # TODO: only for q = 1
    μ₁, μ₂ = N1.μ, N2.μ
    Σ₁, Σ₂ = N1.Σ, N2.Σ
    @assert length(μ₁) == length(μ₂)
    D = length(μ₁)
    return 0.5 * (
        tr(inv(Σ₂)*Σ₁) +
        transpose(μ₂ - μ₁)*inv(Σ₂)*(μ₂ - μ₁) -
        D +
        log(det(Σ₂)/det(Σ₁))
    ) / log(e.base, ℯ)
end

# TODO: define alias divergence, and allow Divergence types to be passed as first argument.
# Can be used with  L2 divergence, for example.

include("analytical.jl")
include("estimators/estimators.jl")
