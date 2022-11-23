export CrossEntropyEstimator
export entropy_cross

"""
    CrossEntropyEstimator <: InformationEstimator

Subtypes of `CrossEntropyEstimator` estimate the differential cross-entropy
(see [`entropy_cross`](@ref)).

Currently implemented subtypes are:

- [`BulinskiDimitrovCE`](@ref).
"""
abstract type CrossEntropyEstimator end

"""
    entropy_cross([e::Entropy,] est::CrossEntropyEstimator,
        x::AbstractDataset{D}, y::AbstractDataset{D})

Estimate the differential cross-entropy of type `e` between `x` and `y`, using
the provided [`CrossEntropyEstimator`](@ref).

The first argument, the entropy type, is optional; it defaults to `Shannon(; base = 2)`.

## Description

## [`Shannon`](@ref) cross entropy

Following the notation of Bulinski & Dimitrov (2021)[^Bulinski2021],
let ``\\mathbb{P}`` and ``\\mathbb{Q}`` be continuous probability measures
with densities ``p(x)`` and ``q(x)``, ``x \\in \\mathcal{R}^D``,
with respect to the Lebesque measure ``\\mu``. Then, writing ``dx := \\mu(dx)``,
Shannon cross entropy is defined as

```math
C(\\mathbb{P}, \\mathbb{Q}) = - \\int_{\\mathbb{R}^d} p(x) \\log{(q(x))} dx.
```

## Other cross entropies

Other types of cross entropies based on entropies where `type(e) != Shannon` are
also possible to compute, but currently no estimators for such quantities are defined
in this package yet.

See also: [`Entropy`](@ref).
"""
function entropy_cross(e::Entropy, est::CrossEntropyEstimator,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D
end

# TODO: we can also define C(ℙ, ℚ) = D(ℙ || ℚ) + H(ℙ), so cross entropy
# can be estimated using any combination of relative entropy and entropy estimators,
# if the estimators are both defined for the same type of entropy (is this true? or
# only for Shannon entropy?)
# Is it wise to mix estimators in this way? Probably not always, but it is technically
# possible.
include("estimators/BulinskiDimitrov.jl")
