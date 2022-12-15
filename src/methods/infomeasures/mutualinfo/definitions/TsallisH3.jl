using Entropies: RectangularBinEncoding, RectangularBinning

export MITsallis
export TsallisH3

"""
    TsallisH3 <: MutualInformationDefinition
    TsallisH3()

A directive to compute the discrete Tsallis mutual *entropy* (Furuichi, 2006)[^TsallisH3].

## Description

Furuichi's Tsallis mutual entropy between variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_q^T(X; Y) = H_q^T(X) - H_q^T(X | Y) = H_q^T(X) + H_q^T(Y) - H_q^T(X, Y),
```

where ``XY \\in \\mathbb{R}^{d_X + d_Y}`` is the joint space, ``H^S(\\cdot)`` and
``H^S(\\cdot, \\cdot)`` are the marginal and joint Shannon entropies, and `q` is the
[`Tsallis`](@ref)-parameter.

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis entropies.
    Journal of Mathematical Physics, 47(2), 023302.

See also: [`mutualinfo`](@ref).
"""
struct TsallisH3 <: MutualInformationDefinition end

# Default behaviour.
estimate(e::Tsallis, est, x, y) = estimate(TsallisH3(), e, est, x, y)

function estimate(def::TsallisH3, measure::MITsallis, est::ProbabilitiesEstimator, x, y)
    e = measure.e
    entropy(e, est, Dataset(x)) +
        entropy(e, est, Dataset(y)) -
        entropy(e, est, Dataset(x, y))
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
function estimate(def::TsallisH3, measure::MITsallis, est::ValueHistogram{RectangularBinning{B}},
        x::AbstractDataset{DX}, y::AbstractDataset{DY}) where {B, DX, DY}
    e = measure.e
    b = est.binning
    XY = Dataset(x, y)
    rb = RectangularBinEncoding(XY, b)
    XY_bins = Dataset([Entropies.encode_as_bin(xyᵢ, rb) for xyᵢ in XY])
    pxy = probabilities(XY_bins)
    px = probabilities(XY_bins[:, 1:DX])
    py = probabilities(XY_bins[:, DX+1:end])
    return entropy(e, px) + entropy(e, py) - entropy(e, pxy)
end


# A different version of Tsallis MI is given in: https://www.mdpi.com/1099-4300/17/8/5382
