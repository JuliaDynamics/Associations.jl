export RelativeEntropyEstimator
export entropy_relative

"""
The supertype for differential relative entropy (see [`entropy_relative`](@ref)) estimators.
Sometimes these are also called Kullback-Leibler (KL) divergence estimators.

Currently implemented subtypes are:

- [`Wang`](@ref).
- [`WangTransformed`](@ref).
"""
abstract type RelativeEntropyEstimator end

"""
    entropy_relative([e::Entropy], est::RelativeEntropyEstimator,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D

Estimate the (differential) relative entropy of type `e` between the `D`-dimensional
datasets `x` and `y` using the provided [`RelativeEntropyEstimator`](@ref).

The first argument, the entropy type, is optional; it defaults to `Shannon(; base = 2)`.

## Description

## [`Shannon`](@ref) relative entropy

We here follow the notation of Bulinski & Dimitrov (2021)[^Bulinski2021].
Let ``\\mathbb{P}`` and ``\\mathbb{Q}`` be continuous probability measures
with densities ``p(x)`` and ``q(x)``, ``x \\in \\mathcal{R}^D``,
with respect to the Lebesque measure ``\\mu``. Then, writing ``dx := \\mu(dx)``, the
Shannon relative entropy is given by

```math
D(\\mathbb{P} || \\mathbb{Q}) =
    \\int_{\\mathbb{R}^D} p(x) \\log{\\left( \\dfrac{p(x)}{q(x)} \\right)} dx,
```

## Other relative entropies

Besides the Shannon relative entropy, one can also define relative entropies based on
other generalized entropies. However, there is no consensus in the literature on
what the definition of such relative entropies are. Therefore, individual
[`RelativeEntropyEstimator`](@ref) specify precisely which quantities they approximate.

See also: [`Entropy`](@ref).

[^Bulinski2021]:
    Bulinski, A., & Dimitrov, D. (2021). Statistical estimation of the Kullback-Leibler
    divergence. Mathematics, 9(5), 544.
"""
function entropy_relative end
# entropy_relative(est::RelativeEntropyEstimator, args...; base = 2, kwargs...) =
#     entropy_relative(Shannon(; base), est, args...; kwargs...)

entropy_relative(e::Entropy, est::RelativeEntropyEstimator, x, y) =
    entropy_relative(e, est, Dataset(x), Dataset(y))
"""
    entropy_relative(::Renyi, p::Probabilities, q::Probabilities)

Estimate the (discrete) relative entropy, or KL divergence, between the pre-computed
probability distributions `p` and `q`, where `p[i]` and `q[i]` is the probability of the
`i`-th outcome in some [outcome_space](@ref) ``\\omega{X}``, defined as

```math
D_{KL}(X || Y) = \\sum_{x \\in \\mathcal{X}} P(x) \\log{\\dfrac{P(x)}{Q(x)}}
```

For this definition to be meaningful `p` and `q` must have *the same* outcome space.

See also: [`probabilities`](@ref).
"""
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

include("analytical.jl")
include("estimators/estimators.jl")
