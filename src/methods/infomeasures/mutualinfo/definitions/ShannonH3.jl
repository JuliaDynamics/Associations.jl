export ShannonH3
"""
    ShannonH3 <: MutualInformationDefinition

A directive to compute the Shannon mutual information ``I^S(X; Y)`` ([`mutualinfo`](@ref))
using a sum of three entropy terms, defined as follows

- Continuous case: ``I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y)``
- Discrete case: ``I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y),

Here, ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint discrete
Shannon entropies, and ``h^S(\\cdot)`` and ``h^S(\\cdot, \\cdot)`` are the corresponding
differential entropies.
"""
struct ShannonH3 <: MutualInformationDefinition end

function estimate(def::ShannonH3, measure::MIShannon, est::EntropyEstimator, x, y)
    e = measure.e
    X = Dataset(x)
    Y = Dataset(y)
    XY = Dataset(X, Y)
    return entropy(e, est, X) + entropy(e, est, Y) - entropy(e, est, XY)
end

function estimate(def::ShannonH3, measure::MIShannon, est::ProbabilitiesEstimator, x, y)
    e = measure.e
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
function estimate(def::ShannonH3, measure::MIShannon, est::ValueHistogram{RectangularBinning{B}},
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
