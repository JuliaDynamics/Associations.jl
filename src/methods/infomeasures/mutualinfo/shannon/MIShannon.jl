export MIShannon
export ShannonH3
"""
    ShannonH3 <: Definition

`ShannonH3` is a directive used in combination with a [`ProbabilitiesEstimator`](@ref) to
compute the discrete Shannon mutual information to base `e.b` between
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

```math
I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y),
```
where ``XY \\in \\mathbb{R}^{d_X + d_Y}`` is the joint space, and ``H^S(\\cdot)`` and
``H^S(\\cdot, \\cdot)`` are the marginal and joint Shannon entropies. Unless otherwise
specified, estimation is done in the most naive way possible: compute probability mass
functions separately for each space `X`, `Y` and `XY`, then plug these probabilites into the
respective entropy formulas.

## Supported estimators

- **[`ValueHistogram`](@ref)**. Bin visitation frequencies are counted in the joint space
    `XY`, then marginal probabilities are obtained from the joint bin visits.
- **[`NaiveKernel``](@ref)**.

See also: [`mutualinfo`](@ref).
"""
struct ShannonH3 <: Definition end

"""
    MIShannon <: MutualInformation
    MIShannon(e = Shannon(; base), definition::Definition = ShannonH3())

A directive to compute the discrete Shannon mutual information using the provided
definition.

## Supported definitions

- [`ShannonH3`](@ref).

## Usage

```julia
using CausalityTools
x, y = Dataset(rand(1000)), Dataset(rand(1000))
b = FixedRectangularBinning(0, 1, 5)
estimate(MIShannon(), ValueHistogram(b), x, y)
```
See also: [`mutualinfo`](@ref).
"""
struct MIShannon{D <: Definition, E <: Renyi} <: MutualInformation
    e::E
    definition::D
    function MIShannon(; base = 2,
            definition::D = ShannonH3()) where {D}
            e = Shannon(; base)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(def::MIShannon{ShannonH3}, est::ProbabilitiesEstimator, x, y)
    e = def.e
    q = e.q
    e.q ≈ 1 || error("$(typeof(est)) estimator not defined for Renyi entropy with q=$(q)")
    X = Dataset(x)
    Y = Dataset(y)
    XY = Dataset(X, Y)
    return entropy(e, est, X) + entropy(e, est, Y) - entropy(e, est, XY)
end

###########################################################################################
# Specialized dispatch.
# ---------------------
# For some information measures and/or methods, specialized dispatch might be necessary
# for estimation to make sense or work to begin with. Define these methods below.
###########################################################################################
function estimate(def::MIShannon{ShannonH3}, est::ValueHistogram{RectangularBinning{B}},
        x::AbstractDataset{DX}, y::AbstractDataset{DY}) where {B, DX, DY}
    e = def.e
    q = e.q
    e.q ≈ 1 || error("$(typeof(est)) estimator not defined for Renyi entropy with q=$(q)")
    b = est.binning

    # Encode joint first, then just sample marginals. This ensures we're using the same
    # bounds for the histogram. For `FixedRectangularBinning` this is not an issue,
    # because bounds are always the same.
    XY = Dataset(x, y)
    rb = RectangularBinEncoding(XY, b)
    XY_bins = Dataset([Entropies.encode_as_bin(xyᵢ, rb) for xyᵢ in XY])
    pxy = probabilities(XY_bins)
    px = probabilities(XY_bins[:, 1:DX])
    py = probabilities(XY_bins[:, DX+1:end])
    return entropy(e, px) + entropy(e, py) - entropy(e, pxy)
end
