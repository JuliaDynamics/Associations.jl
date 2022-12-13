using Entropies: RectangularBinEncoding, RectangularBinning

export MITsallis
export TsallisH3

"""
    TsallisH3 <: Algorithm
    TsallisH3()

`TsallisH3` is a directive to compute the discrete Tsallis mutual *entropy*
introduced by Furuichi (2006)[^TsallisH3].

## Description

Furuichi's Tsallis mutual entropy between variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_q^T(X; Y) = H_q^T(X) - H_q^T(X | Y) = H_q^T(X) + H_q^T(Y) - H_q^T(X, Y),
```

where ``XY \\in \\mathbb{R}^{d_X + d_Y}`` is the joint space, ``H^S(\\cdot)`` and
``H^S(\\cdot, \\cdot)`` are the marginal and joint Shannon entropies, and `q` is the
[`Tsallis`](@ref)-parameter. This is a kind of mutual information
for which certain properties of Shannon mutual information hasn't been proven (see
for reasoning on naming the method). Unless otherwise
    specified, estimation is done in the most naive way possible: compute probability mass
functions separately for each space `X`, `Y` and `XY`, then plug these probabilites into the
respective entropy formulas.

## Supported estimators

- **[`ValueHistogram`](@ref)**. Bin visitation frequencies are counted in the joint space
    `XY`, then marginal probabilities are obtained from the joint bin visits.
- **[`NaiveKernel``](@ref)**.

[^TsallisH3]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis entropies.
    Journal of Mathematical Physics, 47(2), 023302.

See also: [`mutualinfo`](@ref).
"""
struct TsallisH3 <: Definition end

"""
    MITsallis <: MutualInformation
    MITsallis(; base = 2, q = 1.5, k = 0, definition::Definition = TsallisH3())

A directive to compute the discrete Tsallis mutual information using the provided
`definition`.

## Supported definitions

- [`TsallisH3`](@ref).

## Usage

```julia
using CausalityTools
x, y = Dataset(rand(1000, 3)), Dataset(rand(1000, 2))
b = FixedRectangularBinning(0, 1, 5)
est = ValueHistogram(b)
def = MIShannon()
estimate(def, est, x, y)
```
"""
Base.@kwdef struct MITsallis{D <: Definition, E <: Tsallis} <: MutualInformation
    e::E = Tsallis(; q = 1.5, base = 2)
    definition::D = TsallisH3()

    function MITsallis(; base = 2, q = 1.5, k = 0, definition::D = TsallisH3()) where {D}
        e = Tsallis(; q, k, base)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(def::MITsallis{TsallisH3}, est::P, x, y) where P <: ProbabilitiesEstimator
    e = def.e
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
function estimate(def::MITsallis{TsallisH3}, est::ValueHistogram{RectangularBinning{B}},
        x::AbstractDataset{DX}, y::AbstractDataset{DY}) where {B, DX, DY}
    e = def.e
    b = est.binning
    XY = Dataset(x, y)
    rb = RectangularBinEncoding(XY, b)
    XY_bins = Dataset([Entropies.encode_as_bin(xyᵢ, rb) for xyᵢ in XY])
    pxy = probabilities(XY_bins)
    px = probabilities(XY_bins[:, 1:DX])
    py = probabilities(XY_bins[:, DX+1:end])
    return entropy(e, px) + entropy(e, py) - entropy(e, pxy)
end

# """
#     MutualInfoTsallis <: MutualInfoDefinition
#     MutualInfoTsallis()

# `MutualInfoTsallis` is a directive to compute the discrete Tsallis mutual
# information as defined in Tsallis (1998)[^Tsallis1998].

# ## Description

# Assume ``(X, Y)`` is a pair of random variables over the sample space
# ``\\Omega = \\mathcal{X} \\times \\mathcal{Y}``. Let ``P_{XY}(x, y) : \\Omega \\to [0, 1]``
# be the joint distribution of `X` and `Y`, and let ``P_X(x)`` and ``P_Y(y)`` be the marginal
# distributions.

# ```math
# \\begin{align}
# I_{q}^T(X; Y)
# &= D_{q}^T(p(x, y), p(x)p(y)) \\\\
# &= \\dfrac{1}{1 - q}
#     \\left( 1 -
#         \\sum_{x\\in \\mathcal{X}, y\\in \\mathcal{Y}} \\dfrac{P_{XY}(x, y)^q}{P_X(x)^q P_Y(y)^q}
#     \\right)
# \\end{align}
# ```

# [^Tsallis1998]:
#     Tsallis, C. (1998). Generalized entropy-based criterion for consistent testing.
#     Physical Review E, 58(2), 1442.
# """
# struct MutualInfoTsallisTsallis1998 <: MutualInfoDefinition end

# function mutualinfo(e::Tsallis, px, py, pxy) where {P}
#     @assert e.q > 1
#     entropy(e, est.est, Dataset(x)) +
#         entropy(e, est.est, Dataset(y)) -
#         entropy(e, est.est, Dataset(x, y))
# end
