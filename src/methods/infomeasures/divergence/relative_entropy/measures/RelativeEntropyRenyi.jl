
export RelativeEntropyRenyi

"""
    RelativeEntropyRenyi <: Divergence
    RelativeEntropyRenyi(e::Entropy = Shannon())
    RelativeEntropyRenyi(; base::Real)

`RelativeEntropyRenyi` is a directive to compute the discrete Rényi relative entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

- [`RenyiDivergence`](@ref).

See also: [`divergence`](@ref).
"""
struct RelativeEntropyRenyi{E <: Entropy} <: Divergence
    e::E
    function RelativeEntropyRenyi(; base = 2, q = 1.5)
            e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
end

estimate(measure::RelativeEntropyRenyi, est, x, y) =
    estimate(RenyiDivergence(), measure::RelativeEntropyRenyi, est, x, y)

function estimate(def::RenyiDivergence, measure::RelativeEntropyRenyi,
        est::ProbabilitiesEstimator, x, y)
        @show "here"
    q, base = measure.e.q, measure.e.base
    P = probabilities(est, x)
    Q = probabilities(est, y)
    re = 1 / (q - 1) * log(sum(pᵢ^q * qᵢ^(1 - q) for (pᵢ, qᵢ) in zip(P, Q)))
    return re / log(base, ℯ)
end

###########################################################################################
# Specialized dispatch.
# ---------------------
# For some information measures and/or methods, specialized dispatch might be necessary,
# or more efficient. Define these methods below.
###########################################################################################

# Ensures that we're using the same bins for both the joint space and the marginals
# when using `RectangularBinning` (this is automatically the case for
# `FixedRectangularBinning`).
function estimate(def::RenyiDivergence, measure::RelativeEntropyRenyi,
        est::ValueHistogram,
        x::AbstractDataset{DX},
        y::AbstractDataset{DY}) where {DX, DY}
    q, base = measure.e.q, measure.e.base
    b = est.binning
    XY = Dataset(x, y)
    rb = RectangularBinEncoding(XY, b)
    XY_bins = Dataset([Entropies.encode_as_bin(xy, rb) for xy in XY])
    PX = probabilities(XY_bins[:, 1:DX])
    PY = probabilities(XY_bins[:, DX+1:end])
    re = 1 / (q - 1) * log(sum(pᵢ^q * qᵢ^(1 - q) for (pᵢ, qᵢ) in zip(PX, PY)))
    return re / log(base, ℯ)
end


# This function contain analytical expressions for various relative entropies.
# Renyi divergence expressions are from Gil, M. (2011). On Rényi divergence measures for continuous alphabet sources. PhD Thesis.
using Distributions: Beta
using SpecialFunctions: gamma, beta
using Distributions: MvNormal
using LinearAlgebra: tr, inv

function divergence(e::Renyi, x::Beta, y::Beta)
    @assert e.q != 1
    q = e.q # Our q is their α
    αx, αy, βx, βy = x.α, y.α, x.β, y.β
    a = q*αx + (1 - q)*αy
    b = q*βx + (1 - q)*βy

    re = log(beta(αy, βy) / beta(αx, βx)) +
        (1 / (q - 1)) * log(beta(a, b) / beta(αx, βx))
    return re / log(e.base, ℯ)
end

function divergence(e::Renyi, N1::MvNormal, N2::MvNormal)
    @assert e.q != 0
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

# Analytical KL divergence for 1D normals:
# https://ieeexplore.ieee.org/abstract/document/6832827?casa_token=ZhfFH5_G6XgAAAAA:RzQMg0Zjn-CwtOWw4N-jeum3bWzP7tRioSSFAb76fZX58JmXDBW7mSqjbxvr73NDa9fplUSIGw
