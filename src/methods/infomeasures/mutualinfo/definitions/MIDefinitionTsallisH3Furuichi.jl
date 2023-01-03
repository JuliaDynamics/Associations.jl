export MIDefinitionTsallisH3Furuichi

"""
    MIDefinitionTsallisH3Furuichi <: MutualInformationDefinition
    MIDefinitionTsallisH3Furuichi()

A directive to compute the discrete Tsallis mutual *entropy*
(Furuichi, 2006)[^Furuichi2006].

Used in combination with [`MITsallis`](@ref) and [`mutualinfo`](@ref).

## Description

Furuichi's Tsallis mutual entropy between variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_q^T(X; Y) = H_q^T(X) - H_q^T(X | Y) = H_q^T(X) + H_q^T(Y) - H_q^T(X, Y),
```

where ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint Shannon
entropies, and `q` is the [`Tsallis`](@ref)-parameter.

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis entropies.
    Journal of Mathematical Physics, 47(2), 023302.

See also: [`mutualinfo`](@ref).
"""
struct MIDefinitionTsallisH3Furuichi <: MIH3 end
