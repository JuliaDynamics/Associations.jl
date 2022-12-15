"""
    MIH3Shannon <: MutualInformationDefinition

A directive to compute Shannon mutual information as a sum of Shannon entropies.

## Definition (continuous)

```math
I^S(X; Y) = h^S(X) + h^S(Y) - h^S(X, Y),
```

where ``h^S(\\cdot)`` and ``h^S(\\cdot, \\cdot)`` are the marginal and joint differential
Shannon entropies.

## Definition (discrete)

```math
I^S(X; Y) = H^S(X) + H^S(Y) - H_S(X, Y)
```
``h^S(\\cdot)`` and ``h^S(\\cdot, \\cdot)`` are the marginal and joint discrete
Shannon entropies.
"""
struct MIH3Shannon <: MutualInformationDefinition end

"""
    MIH3Tsallis <: MutualInformationDefinition
    MIH3Tsallis()

`MIH3Tsallis` is a directive to compute the discrete Tsallis mutual *entropy*
introduced by Furuichi (2006)[^Furuichi2006].

## Description (discrete)

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

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis entropies.
    Journal of Mathematical Physics, 47(2), 023302.

See also: [`mutualinfo`](@ref).
"""
struct MIH3Tsallis <: MutualInformationDefinition end
