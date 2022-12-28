export ConditionalEntropyTsallis
export Tsallis2H

"""
    Tsallis2H <: ConditionalEntropyDefinition

`Tsallis2H` is a directive used in [`entropy_conditional`](@ref) with a
[`ProbabilityEstimator`](@ref)` to compute the discrete Tsallis conditional entropy
using the formula

```math
H_q^T(Y | X)
= - \\sum_{x, y} p(x, y)^q  \\log_q p(x|y)
= H_q^T(X, Y) - H_q^T(X)
```

where ``\\log_q`` is the `q`-logarithm (see Furuichi, 2006), and `q != 1`.

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis
    entropies. Journal of Mathematical Physics, 47(2), 023302.
"""
struct Tsallis2H <: ConditionalEntropyDefinition end

"""
    ConditionalEntropyTsallis <: ConditionalDifferentialEntropyEstimator
    ConditionalEntropyTsallis(; base = 2, q = 1.5, k = 1,
        definition::Definition = TsallisH2())

`ConditionalEntropyTsallis` is a generic plug-in estimator for discrete conditional
Tsallis entropy.

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

- [`Tsallis2H`](@ref).

## Supported estimators

- [`ValueHistogram`](@ref).

## Usage

```julia
using CausalityTools
x, y = Dataset(rand(1000)), Dataset(rand(1000))
b = FixedRectangularBinning(0, 1, 5)
measure = ConditionalEntropyTsallis(; base = 2, q = 1.7)
entropy_conditional(measure, ValueHistogram(b), x, y)
```
See also: [`mutualinfo`](@ref).
"""
struct ConditionalEntropyTsallis{D <: Definition, E <: Tsallis} <: ConditionalDifferentialEntropyEstimator
    e::E
    definition::D
    function ConditionalEntropyTsallis(; base = 2, q = 1.5, k = 1,
            definition::D = Tsallis2H()) where {D}
            e = Tsallis(; base, q, k)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(measure::ConditionalEntropyTsallis{<:Tsallis2H}, est,
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
function estimate(def::ConditionalEntropyTsallis{Tsallis2H}, est::ValueHistogram{RectangularBinning{B}},
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
