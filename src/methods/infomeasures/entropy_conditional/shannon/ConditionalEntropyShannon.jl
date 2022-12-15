export ConditionalEntropyShannon
export Shannon2H

"""
    Shannon2H <: ConditionalEntropyDefinition

`Shannon2H` is a directive used in [`entropy_conditional`](@ref) with a
[`ProbabilityEstimator`](@ref) to compute the discrete Shannon conditional entropy
using the formula ``H(Y | X) = H(X,Y) - H(X)``.

## Supported estimators

Any [`ProbabilityEstimator`](@ref) that works with multivariate data can in principle
be used. However, not all estimators might be meaningful to use in this context. You should
ensure that the estimator uses the same discretization scheme for both `H(X,Y)` and `H(X)`.
We provide specialized implementations for:

- **[`ValueHistogram`](@ref)**. Bin visitation frequencies are counted in the joint space
    `XY`, then marginal probabilities are obtained from the joint bin visits.
"""
struct Shannon2H <: ConditionalEntropyDefinition end

"""
    ConditionalEntropyShannon <: ConditionalEntropyEstimator
    ConditionalEntropyShannon(; base = 2, definition::Definition = ShannonH2())

`ConditionalEntropyShannon` is a generic plug-in estimator for discrete conditional
Shannon entropy ``H(Y | X)``.

It computes the discrete conditional entropy to the given `base` by first approximating
probabilities using `est`, and plugging them into the formula given by `definition`.
With default settings we simply plug in probabilities to ``H(Y | X) = H(X,Y) - H(X)``.

## Supported definitions

- [`Shannon2H`](@ref).

## Usage

```julia
using CausalityTools
x, y = Dataset(rand(1000)), Dataset(rand(1000))
b = FixedRectangularBinning(0, 1, 5)
measure = ConditionalEntropyShannon(; base = 2)
entropy_conditional(measure, ValueHistogram(b), x, y)
```
See also: [`mutualinfo`](@ref).
"""
struct ConditionalEntropyShannon{D <: Definition, E <: Renyi} <: ConditionalEntropyEstimator
    e::E
    definition::D
    function ConditionalEntropyShannon(; base = 2,
            definition::D = Shannon2H()) where {D}
            e = Shannon(; base)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(measure::ConditionalEntropyShannon{<:Shannon2H}, est,
        x::AbstractDataset, y::AbstractDataset)
    XY = Dataset(x, y)
    X = Dataset(x)
    return entropy(measure.e, est, XY) - entropy(measure.e, est, X)
end

###########################################################################################
# Specialized dispatch.
# ---------------------
# For some information measures and/or methods, specialized dispatch might be necessary
# for estimation to make sense or work to begin with. Define these methods below.
###########################################################################################
function estimate(def::ConditionalEntropyShannon{Shannon2H}, est::ValueHistogram{RectangularBinning{B}},
        x::AbstractDataset{DX}, y::AbstractDataset{DY}) where {B, DX, DY}
    e = def.e
    b = est.binning

    # Encode joint first, then just sample marginals. This ensures we're using the same
    # bounds for the histogram. For `FixedRectangularBinning` this is not an issue,
    # because bounds are always the same.
    XY = Dataset(x, y)
    rb = RectangularBinEncoding(XY, b)
    XY_bins = Dataset([Entropies.encode_as_bin(xyᵢ, rb) for xyᵢ in XY])
    pxy = probabilities(XY_bins)
    px = probabilities(XY_bins[:, 1:DX])
    return entropy(e, pxy) - entropy(e, px)
end
