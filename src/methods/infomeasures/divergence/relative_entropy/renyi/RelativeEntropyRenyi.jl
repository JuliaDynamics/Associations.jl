
export RelativeEntropyRenyi
export RenyiDivergence

"""
    RenyiDivergence  <: DivergenceDefinition
    RenyiDivergence

An instruction to compute the Rényi divergence (relative entropy) according to the
original definition in Rényi (1961)[^Rényi1961]'s seminal paper, here stated in terms of
notation from Van Erven et al. (2014)[[^VanErven2014]].

## Description

For a discrete sample space ``\\Omega`` and probability mass functions
``p(x) : \\Omega \\to [0, 1]`` and ``q(x) : \\Omega \\to [0, 1]``,
the Rényi relative entropy (divergence) is, for `q != 1`, given by

```math
D_q(P || Q) = \\dfrac{1}{q - 1} \\log \\sum_{i = 1}^n p_i^q q_i^{1-\\alpha}.
```

[^Rényi1961]:
    Rényi, A. (1961, June). On measures of entropy and information. In Proceedings of the
    fourth Berkeley symposium on mathematical statistics and probability (Vol. 1, No.
    547-561).

[^VanErven2014]:
    Van Erven, T., & Harremos, P. (2014). Rényi divergence and Kullback-Leibler divergence.
    IEEE Transactions on Information Theory, 60(7), 3797-3820.
"""
struct RenyiDivergence <: DivergenceDefinition end

"""
    RelativeEntropyRenyi <: Divergence
    RelativeEntropyRenyi(e::Entropy = Shannon(), definition = RenyiDivergence())
    RelativeEntropyRenyi(; base::Real)

`RelativeEntropyRenyi` is a directive to compute the discrete Rényi relative entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

- [`RenyiDivergence`](@ref).

See also: [`divergence`](@ref).
"""
struct RelativeEntropyRenyi{D <: DivergenceDefinition, E <: Entropy} <: Divergence
    e::E
    definition::D
    function RelativeEntropyRenyi(; base = 2, q = 1.5,
            definition::D = RenyiDivergence()) where {D}
            e = Renyi(; base, q)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(measure::RelativeEntropyRenyi{<:RenyiDivergence},
        est::ProbabilitiesEstimator, x, y)
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
function estimate(measure::RelativeEntropyRenyi{<:RenyiDivergence},
        est::ValueHistogram,
        x::AbstractDataset{DX},
        y::AbstractDataset{DY}) where {DX, DY}
    q, base = measure.e.q, measure.e.base
    b = est.binning
    XY = Dataset(x, y)
    rb = RectangularBinEncoding(XY, b)
    XY_bins = Dataset([Entropies.encode_as_bin(xyᵢ, rb) for xyᵢ in XY])
    PX = probabilities(XY_bins[:, 1:DX])
    PY = probabilities(XY_bins[:, DX+1:end])
    re = 1 / (q - 1) * log(sum(pᵢ^q * qᵢ^(1 - q) for (pᵢ, qᵢ) in zip(PX, PY)))
    return re / log(base, ℯ)
end
