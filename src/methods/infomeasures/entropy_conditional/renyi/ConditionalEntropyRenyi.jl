export ConditionalEntropyRenyi
export Renyi2H

"""
    Renyi2H <: ConditionalEntropyDefinition

`Renyi2H` is a directive used in [`entropy_conditional`](@ref) with a
[`ProbabilityEstimator`](@ref)` to compute the discrete Rényi conditional entropy
using the following formula (Jizba and Arimutsu, 2004)[^Jizba2004])

```math
R_q(Y | X)
= \\dfrac{1}{q - 1} \\log \\dfrac{\\sum_{x, y} P(x, y)^q}{\\sum_x P(X)^q}
= H^R_q(X, Y) - H^R_q(X).
```

where `q != 1` and ``= H^R_q(X, Y)`` and ``H^R_q(X)`` are joint and marginal Rényi
entropies.

[^Jizba2004]:
    P. Jizba and T. Arimitsu, “The world according to Rényi: Thermodynamics of
    multifractal systems,” Ann. Phys., vol. 312, pp. 17–59, 2004.
"""
struct Renyi2H <: ConditionalEntropyDefinition end

"""
    ConditionalEntropyRenyi <: ConditionalDifferentialEntropyEstimator
    ConditionalEntropyRenyi(; base = 2, q = 1.5,
        definition::Definition = RenyiH2())

`ConditionalEntropyRenyi` is a generic plug-in estimator for discrete conditional
Renyi entropy.

It computes the discrete conditional entropy to the given `base` by first approximating
probabilities using `est`, and plugging them into the formula given by `definition`.
With the default settings, this is a plug-in estimator for

```math
H_q^T(Y | X)
= - \\sum_{x, y} p(x, y)^q  \\log_q p(x|y)
= H_q^T(X, Y) - H_q^T(X)
```

where ``\\log_q`` is the `q`-logarithm (see Furuichi, 2006), and `q != 1`.

## Supported definitions

- [`Renyi2H`](@ref).

## Supported estimators

- [`ValueHistogram`](@ref).

## Usage

```julia
using CausalityTools
x, y = Dataset(rand(1000)), Dataset(rand(1000))
b = FixedRectangularBinning(0, 1, 5)
measure = ConditionalEntropyRenyi(; base = 2, q = 1.7)
entropy_conditional(measure, ValueHistogram(b), x, y)
```
See also: [`mutualinfo`](@ref).
"""
struct ConditionalEntropyRenyi{D <: Definition, E <: Renyi} <: ConditionalDifferentialEntropyEstimator
    e::E
    definition::D
    function ConditionalEntropyRenyi(; base = 2, q = 1.5,
            definition::D = Renyi2H()) where {D}
            e = Renyi(; base, q)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(measure::ConditionalEntropyRenyi{<:Renyi2H}, est,
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
function estimate(def::ConditionalEntropyRenyi{Renyi2H}, est::ValueHistogram{RectangularBinning{B}},
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
